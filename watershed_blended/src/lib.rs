//! Blended-priority metric watershed for non-orthogonal crystal grids.
//!
//! Algorithm: `priority = dgrid_saddle * SCALE - accumulated_path_length`
//!   - dgrid_saddle  : min-propagation (places boundary at physical constriction)
//!   - accumulated_path: sum of physical edge lengths (equalizes flood rate in Å)
//!   - SCALE (default 200): 0.05 Å dgrid resolution per Å of path

use std::collections::BinaryHeap;

use ndarray::Array3;
use numpy::{PyArray3, PyReadonlyArray2, PyReadonlyArray3};
use ordered_float::NotNan;
use pyo3::prelude::*;
// 26-connected neighbor offsets — same order as mgr_rust
// ---------------------------------------------------------------------------
const NEIGHBOR_OFFSETS_26: [(isize, isize, isize); 26] = [
    (-1, -1, -1), (-1, -1, 0), (-1, -1, 1),
    (-1,  0, -1), (-1,  0, 0), (-1,  0, 1),
    (-1,  1, -1), (-1,  1, 0), (-1,  1, 1),
    ( 0, -1, -1), ( 0, -1, 0), ( 0, -1, 1),
    ( 0,  0, -1),                 ( 0,  0, 1),
    ( 0,  1, -1), ( 0,  1, 0), ( 0,  1, 1),
    ( 1, -1, -1), ( 1, -1, 0), ( 1, -1, 1),
    ( 1,  0, -1), ( 1,  0, 0), ( 1,  0, 1),
    ( 1,  1, -1), ( 1,  1, 0), ( 1,  1, 1),
];

// ---------------------------------------------------------------------------
// Compute 26 physical edge lengths from cell matrix + grid shape
// ---------------------------------------------------------------------------
/// Compute 26 physical edge lengths (Å) for the voxel grid.
///
/// `cell_matrix` is ASE `atoms.cell.array` (rows = lattice vectors).
/// `grid_shape` is `(nx, ny, nz)`.
///
/// Step matrix: `step = cell.T @ diag(1/nx, 1/ny, 1/nz)`
/// Edge length for offset `(di,dj,dk)`: `|| step @ [di,dj,dk] ||`
#[pyfunction]
#[pyo3(signature = (cell_matrix, grid_shape))]
fn compute_edge_lengths(
    cell_matrix: PyReadonlyArray2<'_, f64>,
    grid_shape: (usize, usize, usize),
) -> PyResult<Vec<f32>> {
    let cell_arr = cell_matrix.as_array();
    if cell_arr.dim() != (3, 3) {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "cell_matrix must be 3x3",
        ));
    }
    let (nx, ny, nz) = grid_shape;
    if nx == 0 || ny == 0 || nz == 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "grid_shape dimensions must be positive",
        ));
    }

    // Build step matrix: step[i][j] = cell[j][i] / n[j]  (cell^T @ diag(1/n))
    let inv_n = [1.0 / nx as f64, 1.0 / ny as f64, 1.0 / nz as f64];
    let mut step = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            step[i][j] = cell_arr[[j, i]] * inv_n[j]; // cell^T[i][j] * inv_n[j]
        }
    }

    let mut edge_lengths = Vec::with_capacity(26);
    for &(dx, dy, dz) in &NEIGHBOR_OFFSETS_26 {
        let n = [dx as f64, dy as f64, dz as f64];
        // physical displacement = step @ n
        let mut disp = [0.0f64; 3];
        for i in 0..3 {
            for j in 0..3 {
                disp[i] += step[i][j] * n[j];
            }
        }
        let length = (disp[0] * disp[0] + disp[1] * disp[1] + disp[2] * disp[2]).sqrt();
        edge_lengths.push(length as f32);
    }
    Ok(edge_lengths)
}

// ---------------------------------------------------------------------------
// Blended-priority watershed
// ---------------------------------------------------------------------------
/// Blended-priority metric watershed.
///
/// Priority = dgrid_saddle * scale - accumulated_path_length
///
/// - dgrid_saddle  : min-propagation of dgrid values (saddle / constriction detection)
/// - accumulated_path: sum of physical edge lengths from seed (metric correction)
/// - scale          : trades off dgrid resolution vs path correction (default 200)
///
/// Returns uint32 label array with same shape as `image`.
#[pyfunction]
#[pyo3(signature = (image, seeds, edge_lengths, scale=200.0, mask=None))]
fn seeded_watershed_det_blended(
    py: Python<'_>,
    image: PyReadonlyArray3<'_, f32>,
    seeds: PyReadonlyArray3<'_, u32>,
    edge_lengths: Vec<f32>,
    scale: f32,
    mask: Option<PyReadonlyArray3<'_, u8>>,
) -> PyResult<Py<PyArray3<u32>>> {
    let img = image.as_array();
    let seeds_arr = seeds.as_array();
    if img.dim() != seeds_arr.dim() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "image and seeds must have identical shapes",
        ));
    }
    if edge_lengths.len() != 26 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "edge_lengths must have exactly 26 elements (one per 26-NN offset)",
        ));
    }
    if let Some(m) = &mask {
        if m.as_array().dim() != img.dim() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "mask shape must match image",
            ));
        }
    }

    let (nx, ny, nz) = img.dim();
    let total = nx * ny * nz;

    // Flattened arrays for cache-friendly inner loop
    let img_flat: Vec<f32> = img.iter().copied().collect();
    let mask_flat: Option<Vec<u8>> = mask.as_ref().map(|m| {
        m.as_array()
            .iter()
            .map(|&v| if v != 0 { 1u8 } else { 0u8 })
            .collect()
    });

    // Per-voxel state
    let mut labels: Vec<u32> = vec![0; total];
    let mut finalized: Vec<bool> = vec![false; total];
    let mut dgrid_saddle: Vec<f32> = vec![f32::NEG_INFINITY; total];
    let mut path_len: Vec<f32> = vec![0.0; total];
    let mut best_prio: Vec<f32> = vec![f32::NEG_INFINITY; total];

    let to_flat = |i: usize, j: usize, k: usize| -> usize { i * ny * nz + j * nz + k };

    // Max-heap: highest blended priority = closest to seed in (dgrid, path) space
    let mut heap: BinaryHeap<(NotNan<f32>, usize)> =
        BinaryHeap::with_capacity(std::cmp::max(1, total / 8));

    // Seed initialization
    for ((i, j, k), &label) in seeds_arr.indexed_iter() {
        if label == 0 {
            continue;
        }
        let idx = to_flat(i, j, k);
        if let Some(m_flat) = &mask_flat {
            if m_flat[idx] == 0 {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "Seed label found in location where mask==0",
                ));
            }
        }
        labels[idx] = label;
        dgrid_saddle[idx] = img_flat[idx];
        path_len[idx] = 0.0;
        let prio = img_flat[idx] * scale;
        best_prio[idx] = prio;
        heap.push((NotNan::new(prio).unwrap(), idx));
    }

    // Flood
    while let Some((_, idx)) = heap.pop() {
        if finalized[idx] {
            continue;
        }
        finalized[idx] = true;
        // label was set when this entry was pushed

        let i = idx / (ny * nz);
        let j = (idx / nz) % ny;
        let k = idx % nz;

        for (off_idx, &(dx, dy, dz)) in NEIGHBOR_OFFSETS_26.iter().enumerate() {
            let ni = i as isize + dx;
            let nj = j as isize + dy;
            let nk = k as isize + dz;
            if ni < 0 || nj < 0 || nk < 0 {
                continue;
            }
            let (ni_u, nj_u, nk_u) = (ni as usize, nj as usize, nk as usize);
            if ni_u >= nx || nj_u >= ny || nk_u >= nz {
                continue;
            }
            let nidx = to_flat(ni_u, nj_u, nk_u);

            if finalized[nidx] || labels[nidx] != 0 {
                continue;
            }
            if let Some(m_flat) = &mask_flat {
                if m_flat[nidx] == 0 {
                    continue;
                }
            }

            // Min-propagation: dgrid saddle
            let cand_dgrid = img_flat[nidx].min(dgrid_saddle[idx]);
            // Sum-propagation: physical path length
            let cand_path = path_len[idx] + edge_lengths[off_idx];
            // Blended priority
            let cand_prio = cand_dgrid * scale - cand_path;

            if cand_prio > best_prio[nidx] {
                best_prio[nidx] = cand_prio;
                dgrid_saddle[nidx] = cand_dgrid;
                path_len[nidx] = cand_path;
                labels[nidx] = labels[idx];
                heap.push((NotNan::new(cand_prio).unwrap(), nidx));
            }
        }
    }

    let out_arr = Array3::from_shape_vec((nx, ny, nz), labels).expect("shape mismatch");
    Ok(PyArray3::from_owned_array(py, out_arr).to_owned())
}

// ---------------------------------------------------------------------------
// Standard unweighted watershed (for A/B comparison)
// ---------------------------------------------------------------------------
/// Standard unweighted heap watershed (same as mgr_rust::seeded_watershed_det).
/// Provided for side-by-side comparison without importing mgr_rust.
#[pyfunction]
#[pyo3(text_signature = "(image, seeds, mask=None)")]
fn seeded_watershed_det(
    py: Python<'_>,
    image: PyReadonlyArray3<'_, f32>,
    seeds: PyReadonlyArray3<'_, u32>,
    mask: Option<PyReadonlyArray3<'_, u8>>,
) -> PyResult<Py<PyArray3<u32>>> {
    let img = image.as_array();
    let seeds_arr = seeds.as_array();
    if img.dim() != seeds_arr.dim() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "image and seeds must have identical shapes",
        ));
    }
    if let Some(m) = &mask {
        if m.as_array().dim() != img.dim() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "mask shape must match image",
            ));
        }
    }

    let (nx, ny, nz) = img.dim();
    let total = nx * ny * nz;

    let img_flat: Vec<f32> = img.iter().copied().collect();
    let mask_flat: Option<Vec<u8>> = mask.as_ref().map(|m| {
        m.as_array()
            .iter()
            .map(|&v| if v != 0 { 1u8 } else { 0u8 })
            .collect()
    });

    let mut labels: Vec<u32> = vec![0; total];
    let to_flat = |i: usize, j: usize, k: usize| -> usize { i * ny * nz + j * nz + k };

    let mut heap: BinaryHeap<(NotNan<f32>, usize)> =
        BinaryHeap::with_capacity(std::cmp::max(1, total / 8));

    for ((i, j, k), &label) in seeds_arr.indexed_iter() {
        if label == 0 {
            continue;
        }
        let idx = to_flat(i, j, k);
        if let Some(m_flat) = &mask_flat {
            if m_flat[idx] == 0 {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "Seed in mask==0",
                ));
            }
        }
        labels[idx] = label;
        heap.push((NotNan::new(img_flat[idx]).unwrap(), idx));
    }

    while let Some((_prio, idx)) = heap.pop() {
        let i = idx / (ny * nz);
        let j = (idx / nz) % ny;
        let k = idx % nz;
        let current_label = labels[idx];

        for &(dx, dy, dz) in &NEIGHBOR_OFFSETS_26 {
            let ni = i as isize + dx;
            let nj = j as isize + dy;
            let nk = k as isize + dz;
            if ni < 0 || nj < 0 || nk < 0 {
                continue;
            }
            let (ni_u, nj_u, nk_u) = (ni as usize, nj as usize, nk as usize);
            if ni_u >= nx || nj_u >= ny || nk_u >= nz {
                continue;
            }
            let nidx = to_flat(ni_u, nj_u, nk_u);
            if let Some(m_flat) = &mask_flat {
                if m_flat[nidx] == 0 {
                    continue;
                }
            }
            if labels[nidx] == 0 {
                labels[nidx] = current_label;
                heap.push((NotNan::new(img_flat[nidx]).unwrap(), nidx));
            }
        }
    }

    let out_arr = Array3::from_shape_vec((nx, ny, nz), labels).expect("shape mismatch");
    Ok(PyArray3::from_owned_array(py, out_arr).to_owned())
}

// ---------------------------------------------------------------------------
// Module
// ---------------------------------------------------------------------------
#[pymodule]
fn watershed_blended(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_edge_lengths, m)?)?;
    m.add_function(wrap_pyfunction!(seeded_watershed_det_blended, m)?)?;
    m.add_function(wrap_pyfunction!(seeded_watershed_det, m)?)?;
    Ok(())
}
