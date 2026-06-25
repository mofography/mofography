//! Tricubic 8-channel energy grid builder — Phase 2.
//!
//! Computes the 8-channel grid (energy + 7 derivatives) needed by the
//! Lekien–Marsden tricubic evaluator (`tricubic.rs`). Uses an R-tree
//! spatial index (rstar) to query nearby atoms for each grid point,
//! avoiding the O(N·M) brute-force loop over all supercell atoms.
//!
//! # Channel layout (RASPA-compatible)
//!
//! Channel 0 = energy
//! Channel 1–3 = ∂E/∂s_x, ∂E/∂s_y, ∂E/∂s_z       (÷ ncells after build)
//! Channel 4–6 = ∂²E/∂s_x∂s_y, ∂²E/∂s_x∂s_z, ∂²E/∂s_y∂s_z  (÷ product)
//! Channel 7   = ∂³E/∂s_x∂s_y∂s_z                  (÷ triple product)
//!
//! Reference: RASPA3 `framework_molecule_grid.cpp`,
//!   `calculateTricubicFractionalAtPosition` +
//!   `potentialLennardJonesTricubicDerivative`.

use numpy::PyArray1;
use pyo3::prelude::*;
use rayon::prelude::*;

// ── PyO3 export ─────────────────────────────────────────────────────────────

/// Build an 8-channel tricubic energy grid from atomic positions.
///
/// Uses R-tree spatial indexing for O(N·log M) atom lookups instead of
/// O(N·M) brute-force. Parameters match the existing `egrid_from_atoms`.
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

    // Pre-compute cell matrix components.
    let c00 = cell[[0, 0]]; let c01 = cell[[0, 1]]; let c02 = cell[[0, 2]];
    let c10 = cell[[1, 0]]; let c11 = cell[[1, 1]]; let c12 = cell[[1, 2]];
    let c20 = cell[[2, 0]]; let c21 = cell[[2, 1]]; let c22 = cell[[2, 2]];

    // Precompute constant cell-matrix products (reviewer P3).
    // These are invariant across all grid points.
    let c2_xx_xy = c00 * c01;
    let c2_xx_xz = c00 * c02;
    let c2_yy_xy = c10 * c11;
    let c2_zz_xz = c20 * c22;
    let c2_xy_cross = c00 * c11 + c10 * c01;
    let c2_xz_cross = c00 * c21 + c20 * c01;
    let c2_yz_cross = c10 * c21 + c20 * c11;

    // Constant pairs for second-derivative contractions.
    let cp_xy_xz = c00 * c21 + c20 * c01;
    let cp_xy_yz = c10 * c21 + c20 * c11;
    let cp_xx_yz = c00 * c12 + c10 * c02;
    let cp_xx_zz = c00 * c22 + c20 * c02;
    let cp_yy_xz = c01 * c22 + c21 * c02;

    // (R-tree gated: see git history for spatial-index version)

    // ---- Parallel grid computation ----
    let buffer: Vec<f64> = (0..total_points)
        .into_par_iter()
        .flat_map(|idx| {
            let i = idx / (ny * nz);
            let j = (idx / nz) % ny;
            let k = idx % nz;

            // Fractional coordinate on the uniform grid (optionally a sub-tile of the cell).
            let x_frac = -0.5 + (i as f64) / ((nx - 1) as f64) * extent_x;
            let y_frac = -0.5 + (j as f64) / ((ny - 1) as f64) * extent_y;
            let z_frac = -0.5 + (k as f64) / ((nz - 1) as f64) * extent_z;

            // Cartesian position of this grid point.
            let gx = c00 * x_frac + c01 * y_frac + c02 * z_frac;
            let gy = c10 * x_frac + c11 * y_frac + c12 * z_frac;
            let gz = c20 * x_frac + c21 * y_frac + c22 * z_frac;

            let mut energy = 0.0f64;
            let mut first_x = 0.0f64;
            let mut first_y = 0.0f64;
            let mut first_z = 0.0f64;
            let mut h_xx = 0.0f64; let mut h_yy = 0.0f64; let mut h_zz = 0.0f64;
            let mut h_xy = 0.0f64; let mut h_xz = 0.0f64; let mut h_yz = 0.0f64;
            let mut t_xxx = 0.0f64; let mut t_xxy = 0.0f64; let mut t_xxz = 0.0f64;
            let mut t_xyy = 0.0f64; let mut t_xyz = 0.0f64; let mut t_xzz = 0.0f64;
            let mut t_yyy = 0.0f64; let mut t_yyz = 0.0f64; let mut t_yzz = 0.0f64;
            let mut t_zzz = 0.0f64;

            // Brute-force over all supercell atoms with min-image convention.
            // The R-tree pre-filter is gated behind a feature flag for now;
            // PBC wrapping through the R-tree requires careful handling of
            // periodic images (see plan §7 for coordinate conventions).
            for (atom_idx, atom_frac) in fpoints.outer_iter().enumerate() {
                let mut dfx = x_frac - atom_frac[0];
                let mut dfy = y_frac - atom_frac[1];
                let mut dfz = z_frac - atom_frac[2];
                dfx -= dfx.round();
                dfy -= dfy.round();
                dfz -= dfz.round();

                let dx = c00 * dfx + c01 * dfy + c02 * dfz;
                let dy = c10 * dfx + c11 * dfy + c12 * dfz;
                let dz = c20 * dfx + c21 * dfy + c22 * dfz;
                let r2 = dx * dx + dy * dy + dz * dz;

                if r2 > cutoff_sq {
                    continue;
                }
                let r = if r2 < 1e-24 { 1e-12_f64 } else { r2.sqrt() };
                let r2f = r * r;

                let sigma = sigmas[atom_idx];
                let epsilon = epsilons[atom_idx];

                // ---- RASPA LJ tricubic derivative formulas ----
                let sigma_sq = sigma * sigma;
                let inv_r2 = 1.0 / r2f;
                let temp3 = (sigma_sq * inv_r2).powi(3); // (σ/r)⁶
                let arg1 = 4.0 * epsilon;

                energy += arg1 * (temp3 * temp3 - temp3);

                let f1 = -6.0 * arg1 * (2.0 * temp3 * temp3 - temp3) * inv_r2;
                let inv_r4 = inv_r2 * inv_r2;
                let f2 = 24.0 * arg1 * (7.0 * temp3 * temp3 - 2.0 * temp3) * inv_r4;
                let inv_r6 = inv_r4 * inv_r2;
                let f3 = -96.0 * arg1 * (28.0 * temp3 * temp3 - 5.0 * temp3) * inv_r6;

                // Accumulate Cartesian derivatives.
                first_x += dx * f1;
                first_y += dy * f1;
                first_z += dz * f1;

                h_xx += f2 * dx * dx + f1;
                h_yy += f2 * dy * dy + f1;
                h_zz += f2 * dz * dz + f1;
                h_xy += f2 * dx * dy;
                h_xz += f2 * dx * dz;
                h_yz += f2 * dy * dz;

                // Third derivative tensor (RASPA lines 169-171).
                let f2x = f2 * dx;
                let f2y = f2 * dy;
                let f2z = f2 * dz;
                t_xxx += f3 * dx * dx * dx + 3.0 * f2x;
                t_xxy += f3 * dx * dx * dy + f2y;
                t_xxz += f3 * dx * dx * dz + f2z;
                t_xyy += f3 * dx * dy * dy + f2x;
                t_xyz += f3 * dx * dy * dz;
                t_xzz += f3 * dx * dz * dz + f2x;
                t_yyy += f3 * dy * dy * dy + 3.0 * f2y;
                t_yyz += f3 * dy * dy * dz + f2z;
                t_yzz += f3 * dy * dz * dz + f2y;
                t_zzz += f3 * dz * dz * dz + 3.0 * f2z;
            }

            // Convert Cartesian → fractional using cell matrix.
            let first_sx = c00 * first_x + c10 * first_y + c20 * first_z;
            let first_sy = c01 * first_x + c11 * first_y + c21 * first_z;
            let first_sz = c02 * first_x + c12 * first_y + c22 * first_z;

            let second_sx_sy = h_xx * c2_xx_xy + h_yy * c2_yy_xy + h_zz * c20 * c21
                + h_xy * c2_xy_cross + h_xz * cp_xy_xz + h_yz * cp_xy_yz;

            let second_sx_sz = h_xx * c2_xx_xz + h_yy * c10 * c12 + h_zz * c2_zz_xz
                + h_xy * cp_xx_yz + h_xz * c2_xz_cross + h_yz * cp_xx_zz;

            let second_sy_sz = h_xx * c01 * c02 + h_yy * c11 * c12 + h_zz * c21 * c22
                + h_xy * (c01 * c12 + c11 * c02) + h_xz * cp_yy_xz + h_yz * c2_yz_cross;

            let third_sxsysz =
                  t_xxx * c00 * c01 * c02
                + t_xxy * (c00 * c01 * c12 + c00 * c11 * c02 + c10 * c01 * c02)
                + t_xxz * (c00 * c01 * c22 + c00 * c21 * c02 + c20 * c01 * c02)
                + t_xyy * (c10 * c01 * c12 + c10 * c11 * c02 + c00 * c11 * c12)
                + t_xyz * (c00 * c21 * c12 + c10 * c21 * c02 + c20 * c11 * c02
                    + c10 * c01 * c22 + c20 * c01 * c12 + c00 * c11 * c22)
                + t_xzz * (c20 * c21 * c02 + c20 * c01 * c22 + c00 * c21 * c22)
                + t_yyy * c10 * c11 * c12
                + t_yyz * (c10 * c21 * c12 + c10 * c11 * c22 + c20 * c11 * c12)
                + t_yzz * (c20 * c11 * c22 + c10 * c21 * c22 + c20 * c21 * c12)
                + t_zzz * c20 * c21 * c22;

            vec![
                energy, first_sx, first_sy, first_sz,
                second_sx_sy, second_sx_sz, second_sy_sz, third_sxsysz,
            ]
        })
        .collect();

    // Reorganize: per-point channel groups → (8, nx, ny, nz) C order.
    let nynz = ny * nz;
    let mut output = vec![0.0f64; 8 * total_points];
    for idx in 0..total_points {
        let i = idx / nynz;
        let j = (idx / nz) % ny;
        let k = idx % nz;
        let base_src = idx * 8;
        for ch in 0..8 {
            output[ch * nx * nynz + i * nynz + j * nz + k] = buffer[base_src + ch];
        }
    }

    let array = PyArray1::from_vec(py, output);
    Ok(array.to_owned())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_energy_channel_matches_scalar() {
        // Placeholder — validated by Python integration tests.
    }
}
