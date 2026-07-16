//! 8-channel Coulomb potential grid builder (Wolf DSF, φ_K in Kelvin).

use crate::coulomb_wolf::{accumulate_wolf_dsf, WolfConstants};
use crate::egrid_tricubic_common::{
    frac_delta_to_cart, grid_fractional_coords, min_image_frac_delta, reorganize_to_channel_major,
    CartesianDerivatives, CellContraction,
};
use numpy::PyArray1;
use pyo3::prelude::*;
use rayon::prelude::*;

/// Build an 8-channel Coulomb potential grid φ_K(r) from framework atom charges.
///
/// Uses Wolf DSF. Channel 0 = φ_K (Kelvin per unit probe charge +1 e).
/// Multiply by guest charge at interpolation: E_coul = q_guest × φ_K.
#[pyfunction]
#[pyo3(signature = (fpoints, cell, charges, nx, ny, nz, cutoff, alpha=None, extent_x=1.0, extent_y=1.0, extent_z=1.0))]
#[pyo3(text_signature = "(fpoints, cell, charges, nx, ny, nz, cutoff, alpha=None, extent_x=1.0, extent_y=1.0, extent_z=1.0)")]
pub fn coulomb_tricubic_from_atoms(
    py: Python<'_>,
    fpoints: numpy::PyReadonlyArray2<'_, f64>,
    cell: numpy::PyReadonlyArray2<'_, f64>,
    charges: numpy::PyReadonlyArray1<'_, f64>,
    nx: usize,
    ny: usize,
    nz: usize,
    cutoff: f64,
    alpha: Option<f64>,
    extent_x: f64,
    extent_y: f64,
    extent_z: f64,
) -> PyResult<Py<PyArray1<f64>>> {
    if nx < 2 || ny < 2 || nz < 2 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "nx, ny, nz must each be >= 2 for tricubic grids",
        ));
    }

    let fpoints = fpoints.as_array();
    let cell = cell.as_array();
    let charges = charges.as_array();

    if fpoints.nrows() != charges.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "charges length must match fpoints rows",
        ));
    }

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

    let alpha = alpha.unwrap_or(4.0 / cutoff);
    let w = WolfConstants::new(alpha, cutoff);
    let cutoff_sq = cutoff * cutoff;
    let total_points = nx * ny * nz;

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

                let q_j = charges[atom_idx];
                let (mut f1, mut f2, mut f3) = (0.0, 0.0, 0.0);
                accumulate_wolf_dsf(
                    &w,
                    q_j,
                    r,
                    r2,
                    &mut derivs.energy,
                    &mut f1,
                    &mut f2,
                    &mut f3,
                );
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
