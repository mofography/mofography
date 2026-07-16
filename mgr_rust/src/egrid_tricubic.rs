//! Tricubic 8-channel energy grid builder — Phase 2.
//!
//! Computes the 8-channel grid (energy + 7 derivatives) needed by the
//! Lekien–Marsden tricubic evaluator (`tricubic.rs`).
//!
//! # Channel layout (RASPA-compatible)
//!
//! Channel 0 = energy (Kelvin, ε/k_B)
//! Channel 1–3 = ∂E/∂s_x, ∂E/∂s_y, ∂E/∂s_z       (grid-step units)
//! Channel 4–6 = ∂²E/∂s_x∂s_y, ∂²E/∂s_x∂s_z, ∂²E/∂s_y∂s_z  (grid-step units)
//! Channel 7   = ∂³E/∂s_x∂s_y∂s_z                  (grid-step units)

use crate::egrid_tricubic_common::{
    frac_delta_to_cart, grid_fractional_coords, min_image_frac_delta, reorganize_to_channel_major,
    CartesianDerivatives, CellContraction,
};
use numpy::PyArray1;
use pyo3::prelude::*;
use rayon::prelude::*;

/// Build an 8-channel tricubic energy grid from atomic positions.
///
/// `fpoints` should already be supercell-replicated by the Python wrapper.
#[pyfunction]
#[pyo3(signature = (fpoints, cell, epsilons, sigmas, nx, ny, nz, cutoff, extent_x=1.0, extent_y=1.0, extent_z=1.0))]
#[pyo3(text_signature = "(fpoints, cell, epsilons, sigmas, nx, ny, nz, cutoff, extent_x=1.0, extent_y=1.0, extent_z=1.0)")]
pub fn egrid_tricubic_from_atoms(
    py: Python<'_>,
    fpoints: numpy::PyReadonlyArray2<'_, f64>,
    cell: numpy::PyReadonlyArray2<'_, f64>,
    epsilons: numpy::PyReadonlyArray1<'_, f64>,
    sigmas: numpy::PyReadonlyArray1<'_, f64>,
    nx: usize,
    ny: usize,
    nz: usize,
    cutoff: f64,
    extent_x: f64,
    extent_y: f64,
    extent_z: f64,
) -> PyResult<Py<PyArray1<f64>>> {
    let fpoints = fpoints.as_array();
    let cell = cell.as_array();
    let epsilons = epsilons.as_array();
    let sigmas = sigmas.as_array();

    let total_points = nx * ny * nz;
    let cutoff_sq = cutoff * cutoff;

    let c00 = cell[[0, 0]];
    let c01 = cell[[0, 1]];
    let c02 = cell[[0, 2]];
    let c10 = cell[[1, 0]];
    let c11 = cell[[1, 1]];
    let c12 = cell[[1, 2]];
    let c20 = cell[[2, 0]];
    let c21 = cell[[2, 1]];
    let c22 = cell[[2, 2]];
    let cell_contract = CellContraction::from_cell_matrix(
        c00, c01, c02, c10, c11, c12, c20, c21, c22,
    );

    let buffer: Vec<f64> = (0..total_points)
        .into_par_iter()
        .flat_map(|idx| {
            let i = idx / (ny * nz);
            let j = (idx / nz) % ny;
            let k = idx % nz;

            let (x_frac, y_frac, z_frac) =
                grid_fractional_coords(i, j, k, nx, ny, nz, extent_x, extent_y, extent_z);

            let mut derivs = CartesianDerivatives::default();

            for (atom_idx, atom_frac) in fpoints.outer_iter().enumerate() {
                let (dfx, dfy, dfz) = min_image_frac_delta(
                    x_frac,
                    y_frac,
                    z_frac,
                    atom_frac[0],
                    atom_frac[1],
                    atom_frac[2],
                );
                let (dx, dy, dz) = frac_delta_to_cart(&cell_contract, dfx, dfy, dfz);
                let r2 = dx * dx + dy * dy + dz * dz;

                if r2 > cutoff_sq {
                    continue;
                }
                let r = if r2 < 1e-24 { 1e-12_f64 } else { r2.sqrt() };
                let r2f = r * r;

                let sigma = sigmas[atom_idx];
                let epsilon = epsilons[atom_idx];

                let sigma_sq = sigma * sigma;
                let inv_r2 = 1.0 / r2f;
                let temp3 = (sigma_sq * inv_r2).powi(3);
                let arg1 = 4.0 * epsilon;

                derivs.energy += arg1 * (temp3 * temp3 - temp3);

                let f1 = -6.0 * arg1 * (2.0 * temp3 * temp3 - temp3) * inv_r2;
                let inv_r4 = inv_r2 * inv_r2;
                let f2 = 24.0 * arg1 * (7.0 * temp3 * temp3 - 2.0 * temp3) * inv_r4;
                let inv_r6 = inv_r4 * inv_r2;
                let f3 = -96.0 * arg1 * (28.0 * temp3 * temp3 - 5.0 * temp3) * inv_r6;

                derivs.accumulate_radial(dx, dy, dz, f1, f2, f3);
            }

            derivs
                .to_fractional_channels(&cell_contract, nx, ny, nz, extent_x, extent_y, extent_z)
                .to_vec()
        })
        .collect();

    let output = reorganize_to_channel_major(buffer, nx, ny, nz);
    let array = PyArray1::from_vec(py, output);
    Ok(array.to_owned())
}
