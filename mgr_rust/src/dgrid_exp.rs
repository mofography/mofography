//! Optimized `dgrid_from_atoms`: SoA layout, SIMD (orthogonal), 7-image MIC on skewed cells
//! followed by an exact refinement of the voxels the 7-image rule cannot certify.
//!
//! The 7-image set (rounded image plus six face shifts) is not an exact minimum image on
//! skewed cells: the nearest image can need an edge or corner shift such as (0,-1,1). A
//! rounded image is exact whenever the true separation is below half the smallest
//! perpendicular cell width `h_min` (its fractional components then lie inside (-1/2, 1/2)).
//! A voxel whose result `d` satisfies `d + r_max < h_min / 2` is therefore certified exact.
//! It is also certified when only one axis can shift (`d + r_max >= h_i / 2` for a single
//! axis), because the face shifts along that axis are already tested. Every other voxel of a
//! skewed cell is recomputed over {-1,0,1}^3 (or {-2..2}^3 when a sampled check shows the
//! 27 shifts are not enough).

use ndarray::Array3;
use numpy::{PyArray3, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::Py;
use rayon::prelude::*;
use wide::f64x4;

const ORTHO_TOL: f64 = 1e-12;

const OFFSETS_FACE: [[f64; 3]; 6] = [
    [1.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, -1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, -1.0],
];

#[derive(Clone, Copy)]
struct CellMat {
    c00: f64,
    c01: f64,
    c02: f64,
    c10: f64,
    c11: f64,
    c12: f64,
    c20: f64,
    c21: f64,
    c22: f64,
}

struct AtomsSoA {
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    r: Vec<f64>,
}

impl AtomsSoA {
    fn from_views(
        fpoints: &ndarray::ArrayView2<f64>,
        radii: &ndarray::ArrayView1<f64>,
    ) -> Self {
        let n = fpoints.len_of(ndarray::Axis(0));
        let mut x = Vec::with_capacity(n);
        let mut y = Vec::with_capacity(n);
        let mut z = Vec::with_capacity(n);
        let mut r = Vec::with_capacity(n);
        for (idx, atom_frac) in fpoints.outer_iter().enumerate() {
            x.push(atom_frac[0]);
            y.push(atom_frac[1]);
            z.push(atom_frac[2]);
            r.push(radii[idx]);
        }
        Self { x, y, z, r }
    }

    fn len(&self) -> usize {
        self.x.len()
    }
}

#[inline]
fn cell_is_orthogonal(cell: &ndarray::ArrayView2<f64>) -> bool {
    cell[[0, 1]].abs() < ORTHO_TOL
        && cell[[0, 2]].abs() < ORTHO_TOL
        && cell[[1, 0]].abs() < ORTHO_TOL
        && cell[[1, 2]].abs() < ORTHO_TOL
        && cell[[2, 0]].abs() < ORTHO_TOL
        && cell[[2, 1]].abs() < ORTHO_TOL
}

#[inline]
fn d_max_sq_from_cell(cell: &ndarray::ArrayView2<f64>) -> f64 {
    let col_len_sq = |col: usize| {
        let c0 = cell[[0, col]];
        let c1 = cell[[1, col]];
        let c2 = cell[[2, col]];
        c0 * c0 + c1 * c1 + c2 * c2
    };
    let min_len_sq = col_len_sq(0).min(col_len_sq(1)).min(col_len_sq(2));
    0.25 * min_len_sq
}

#[inline]
fn frac_to_cart_sq(dx: f64, dy: f64, dz: f64, cell: CellMat) -> f64 {
    let cart_x = cell.c00 * dx + cell.c01 * dy + cell.c02 * dz;
    let cart_y = cell.c10 * dx + cell.c11 * dy + cell.c12 * dz;
    let cart_z = cell.c20 * dx + cell.c21 * dy + cell.c22 * dz;
    cart_x * cart_x + cart_y * cart_y + cart_z * cart_z
}

#[inline]
fn maybe_update_min_dist(min_dist: &mut f64, sq_dist: f64, r_atom: f64) {
    let threshold = *min_dist + r_atom;
    if threshold > 0.0 && sq_dist >= threshold * threshold {
        return;
    }
    let dist = sq_dist.sqrt() - r_atom;
    if dist < *min_dist {
        *min_dist = dist;
    }
}

#[inline]
fn mic_min_dist_skewed(
    raw_dx: f64,
    raw_dy: f64,
    raw_dz: f64,
    r_atom: f64,
    cell: CellMat,
    d_max_sq: f64,
    min_dist: &mut f64,
) {
    let dx = raw_dx - raw_dx.round();
    let dy = raw_dy - raw_dy.round();
    let dz = raw_dz - raw_dz.round();
    let sq_nint = frac_to_cart_sq(dx, dy, dz, cell);

    if sq_nint < d_max_sq {
        maybe_update_min_dist(min_dist, sq_nint, r_atom);
        return;
    }

    maybe_update_min_dist(min_dist, sq_nint, r_atom);
    for off in OFFSETS_FACE {
        let sq = frac_to_cart_sq(dx + off[0], dy + off[1], dz + off[2], cell);
        maybe_update_min_dist(min_dist, sq, r_atom);
    }
}

#[inline]
fn mic_min_dist_orthogonal_simd4(
    x_frac: f64,
    y_frac: f64,
    z_frac: f64,
    atoms: &AtomsSoA,
    base: usize,
    cell: CellMat,
    min_dist: &mut f64,
) {
    let xf = f64x4::splat(x_frac);
    let yf = f64x4::splat(y_frac);
    let zf = f64x4::splat(z_frac);

    let rdx = xf
        - f64x4::new([
            atoms.x[base],
            atoms.x[base + 1],
            atoms.x[base + 2],
            atoms.x[base + 3],
        ]);
    let rdy = yf
        - f64x4::new([
            atoms.y[base],
            atoms.y[base + 1],
            atoms.y[base + 2],
            atoms.y[base + 3],
        ]);
    let rdz = zf
        - f64x4::new([
            atoms.z[base],
            atoms.z[base + 1],
            atoms.z[base + 2],
            atoms.z[base + 3],
        ]);

    let dx = rdx - rdx.round();
    let dy = rdy - rdy.round();
    let dz = rdz - rdz.round();

    let c00 = f64x4::splat(cell.c00);
    let c01 = f64x4::splat(cell.c01);
    let c02 = f64x4::splat(cell.c02);
    let c10 = f64x4::splat(cell.c10);
    let c11 = f64x4::splat(cell.c11);
    let c12 = f64x4::splat(cell.c12);
    let c20 = f64x4::splat(cell.c20);
    let c21 = f64x4::splat(cell.c21);
    let c22 = f64x4::splat(cell.c22);

    let cart_x = c00 * dx + c01 * dy + c02 * dz;
    let cart_y = c10 * dx + c11 * dy + c12 * dz;
    let cart_z = c20 * dx + c21 * dy + c22 * dz;
    let sq = cart_x * cart_x + cart_y * cart_y + cart_z * cart_z;

    let sq_arr = sq.to_array();
    for lane in 0..4 {
        maybe_update_min_dist(min_dist, sq_arr[lane], atoms.r[base + lane]);
    }
}

#[inline]
fn min_dist_at_voxel(
    x_frac: f64,
    y_frac: f64,
    z_frac: f64,
    atoms: &AtomsSoA,
    cell: CellMat,
    is_orthogonal: bool,
    d_max_sq: f64,
) -> f32 {
    let mut min_dist = f64::MAX;
    let n = atoms.len();
    let mut i = 0;

    if is_orthogonal {
        while i + 4 <= n {
            mic_min_dist_orthogonal_simd4(x_frac, y_frac, z_frac, atoms, i, cell, &mut min_dist);
            i += 4;
        }
    }

    while i < n {
        let raw_dx = x_frac - atoms.x[i];
        let raw_dy = y_frac - atoms.y[i];
        let raw_dz = z_frac - atoms.z[i];
        let r_atom = atoms.r[i];

        if is_orthogonal {
            let dx = raw_dx - raw_dx.round();
            let dy = raw_dy - raw_dy.round();
            let dz = raw_dz - raw_dz.round();
            let sq = frac_to_cart_sq(dx, dy, dz, cell);
            maybe_update_min_dist(&mut min_dist, sq, r_atom);
        } else {
            mic_min_dist_skewed(raw_dx, raw_dy, raw_dz, r_atom, cell, d_max_sq, &mut min_dist);
        }
        i += 1;
    }

    min_dist as f32
}

/// All integer shifts in {-n..n}^3.
fn image_shifts(n: i32) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    for i in -n..=n {
        for j in -n..=n {
            for k in -n..=n {
                out.push([i as f64, j as f64, k as f64]);
            }
        }
    }
    out
}

/// Perpendicular width of the cell along each axis (distance between opposite faces).
fn perp_widths(cell: &ndarray::ArrayView2<f64>) -> [f64; 3] {
    let col = |c: usize| [cell[[0, c]], cell[[1, c]], cell[[2, c]]];
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let norm = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let (a, b, c) = (col(0), col(1), col(2));
    let bc = cross(b, c);
    let vol = (a[0] * bc[0] + a[1] * bc[1] + a[2] * bc[2]).abs();
    [
        vol / norm(bc),
        vol / norm(cross(c, a)),
        vol / norm(cross(a, b)),
    ]
}

/// True when some separation is shorter over {-2..2}^3 than over {-1,0,1}^3, judged on a
/// fixed pseudo-random sample of reduced fractional vectors. Selects the 125-shift refine.
fn needs_wide_shell(cell: CellMat) -> bool {
    let near = image_shifts(1);
    let wide = image_shifts(2);
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    let best = |g: [f64; 3], shifts: &[[f64; 3]]| {
        shifts
            .iter()
            .map(|s| frac_to_cart_sq(g[0] + s[0], g[1] + s[1], g[2] + s[2], cell))
            .fold(f64::MAX, f64::min)
    };
    (0..4096).any(|_| {
        let g = [next(), next(), next()];
        best(g, &wide) < best(g, &near) * (1.0 - 1e-12)
    })
}

/// Exact minimum over the shifts {-range..range}^3, starting from the 7-image result
/// `d_start` (itself the distance to a real image, so a valid upper bound). Along axis `i`
/// an image with shift `s_i` is at least `|g_i + s_i| * widths[i]` away, so shifts that
/// cannot beat the current minimum are never converted to Cartesian.
#[inline]
fn min_dist_at_voxel_images(
    x_frac: f64,
    y_frac: f64,
    z_frac: f64,
    atoms: &AtomsSoA,
    cell: CellMat,
    widths: [f64; 3],
    range: i32,
    d_start: f32,
) -> f32 {
    let mut min_dist = d_start as f64;
    for i in 0..atoms.len() {
        let r_atom = atoms.r[i];
        let reach = min_dist + r_atom;
        if reach <= 0.0 {
            continue;
        }
        let raw = [x_frac - atoms.x[i], y_frac - atoms.y[i], z_frac - atoms.z[i]];
        let g = [raw[0] - raw[0].round(), raw[1] - raw[1].round(), raw[2] - raw[2].round()];
        let mut list = [[0.0f64; 5]; 3];
        let mut len = [0usize; 3];
        for ax in 0..3 {
            for sh in -range..=range {
                let s = sh as f64;
                if (g[ax] + s).abs() * widths[ax] < reach {
                    list[ax][len[ax]] = s;
                    len[ax] += 1;
                }
            }
        }
        for ix in 0..len[0] {
            for iy in 0..len[1] {
                for iz in 0..len[2] {
                    let sq = frac_to_cart_sq(
                        g[0] + list[0][ix],
                        g[1] + list[1][iy],
                        g[2] + list[2][iz],
                        cell,
                    );
                    maybe_update_min_dist(&mut min_dist, sq, r_atom);
                }
            }
        }
    }
    min_dist as f32
}

/// Work done by the skewed-cell refinement, for benchmarking.
#[derive(Clone, Copy, Default)]
pub struct DgridStats {
    pub total_voxels: usize,
    pub refined_voxels: usize,
    pub changed_voxels: usize,
    pub max_change: f32,
    pub main_seconds: f64,
    pub refine_seconds: f64,
    pub wide_shell: bool,
}

fn dgrid_impl(
    fpoints: &ndarray::ArrayView2<f64>,
    cell_view: &ndarray::ArrayView2<f64>,
    radii: &ndarray::ArrayView1<f64>,
    nx: usize,
    ny: usize,
    nz: usize,
    refine: bool,
) -> (Array3<f32>, DgridStats) {
    let atoms = AtomsSoA::from_views(fpoints, radii);
    let total_points = nx * ny * nz;

    let cell_mat = CellMat {
        c00: cell_view[[0, 0]],
        c01: cell_view[[0, 1]],
        c02: cell_view[[0, 2]],
        c10: cell_view[[1, 0]],
        c11: cell_view[[1, 1]],
        c12: cell_view[[1, 2]],
        c20: cell_view[[2, 0]],
        c21: cell_view[[2, 1]],
        c22: cell_view[[2, 2]],
    };

    let is_orthogonal = cell_is_orthogonal(cell_view);
    let d_max_sq = if is_orthogonal {
        0.0
    } else {
        d_max_sq_from_cell(cell_view)
    };

    let frac_of = move |idx: usize| {
        let i = idx / (ny * nz);
        let j = (idx / nz) % ny;
        let k = idx % nz;
        (
            -0.5 + (i as f64) / ((nx - 1) as f64),
            -0.5 + (j as f64) / ((ny - 1) as f64),
            -0.5 + (k as f64) / ((nz - 1) as f64),
        )
    };

    let t_main = std::time::Instant::now();
    let mut buffer: Vec<f32> = (0..total_points)
        .into_par_iter()
        .map(|idx| {
            let (x_frac, y_frac, z_frac) = frac_of(idx);
            min_dist_at_voxel(
                x_frac,
                y_frac,
                z_frac,
                &atoms,
                cell_mat,
                is_orthogonal,
                d_max_sq,
            )
        })
        .collect();

    let mut stats = DgridStats {
        total_voxels: total_points,
        main_seconds: t_main.elapsed().as_secs_f64(),
        ..DgridStats::default()
    };

    if refine && !is_orthogonal && atoms.len() > 0 {
        let t_refine = std::time::Instant::now();
        let r_max = atoms.r.iter().cloned().fold(f64::MIN, f64::max);
        let widths = perp_widths(cell_view);
        // An image shifted along axis i is at least (1 - |g_i|) * widths[i] >= widths[i] / 2
        // away, and the true nearest pair is within d + r_max. A voxel is certified when at
        // most one axis can shift (face shifts, already tested, cover that) and no
        // double shift (>= 1.5 * widths[i]) is in reach.
        let reach_tol = 1.0 - 1e-9;
        let uncertified: Vec<usize> = buffer
            .par_iter()
            .enumerate()
            .filter(|(_, &d)| {
                let reach = (d as f64) + r_max;
                let single = widths.iter().filter(|w| reach >= 0.5 * **w * reach_tol).count();
                single >= 2 || widths.iter().any(|w| reach >= 1.5 * *w * reach_tol)
            })
            .map(|(idx, _)| idx)
            .collect();
        stats.refined_voxels = uncertified.len();
        if !uncertified.is_empty() {
            stats.wide_shell = needs_wide_shell(cell_mat);
            let range = if stats.wide_shell { 2 } else { 1 };
            let buffer_ref = &buffer;
            let refined: Vec<f32> = uncertified
                .par_iter()
                .map(|&idx| {
                    let (x_frac, y_frac, z_frac) = frac_of(idx);
                    min_dist_at_voxel_images(
                        x_frac, y_frac, z_frac, &atoms, cell_mat, widths, range, buffer_ref[idx],
                    )
                })
                .collect();
            for (&idx, &d_new) in uncertified.iter().zip(refined.iter()) {
                let change = buffer[idx] - d_new;
                if change > 1e-6 {
                    stats.changed_voxels += 1;
                    stats.max_change = stats.max_change.max(change);
                }
                buffer[idx] = d_new.min(buffer[idx]);
            }
        }
        stats.refine_seconds = t_refine.elapsed().as_secs_f64();
    }

    let array = Array3::from_shape_vec((nx, ny, nz), buffer)
        .expect("Shape mismatch while reshaping the output distance grid");
    (array, stats)
}

#[pyfunction]
#[pyo3(text_signature = "(fpoints, cell, radii, nx, ny, nz)")]
pub fn dgrid_from_atoms(
    py: Python<'_>,
    fpoints: PyReadonlyArray2<'_, f64>,
    cell: PyReadonlyArray2<'_, f64>,
    radii: PyReadonlyArray1<'_, f64>,
    nx: usize,
    ny: usize,
    nz: usize,
) -> PyResult<Py<PyArray3<f32>>> {
    let (array, _) = dgrid_impl(
        &fpoints.as_array(),
        &cell.as_array(),
        &radii.as_array(),
        nx,
        ny,
        nz,
        true,
    );
    Ok(PyArray3::from_owned_array(py, array).to_owned())
}

/// Same grid as `dgrid_from_atoms`, plus refinement statistics:
/// `(grid, {total_voxels, refined_voxels, changed_voxels, max_change, main_seconds,
/// refine_seconds, wide_shell})`. `refine=False` returns the raw 7-image grid.
#[pyfunction]
#[pyo3(signature = (fpoints, cell, radii, nx, ny, nz, refine=true))]
pub fn dgrid_from_atoms_stats(
    py: Python<'_>,
    fpoints: PyReadonlyArray2<'_, f64>,
    cell: PyReadonlyArray2<'_, f64>,
    radii: PyReadonlyArray1<'_, f64>,
    nx: usize,
    ny: usize,
    nz: usize,
    refine: bool,
) -> PyResult<(Py<PyArray3<f32>>, PyObject)> {
    let (array, st) = dgrid_impl(
        &fpoints.as_array(),
        &cell.as_array(),
        &radii.as_array(),
        nx,
        ny,
        nz,
        refine,
    );
    let d = pyo3::types::PyDict::new(py);
    d.set_item("total_voxels", st.total_voxels)?;
    d.set_item("refined_voxels", st.refined_voxels)?;
    d.set_item("changed_voxels", st.changed_voxels)?;
    d.set_item("max_change", st.max_change)?;
    d.set_item("main_seconds", st.main_seconds)?;
    d.set_item("refine_seconds", st.refine_seconds)?;
    d.set_item("wide_shell", st.wide_shell)?;
    Ok((PyArray3::from_owned_array(py, array).to_owned(), d.into()))
}
