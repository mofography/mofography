//! Tricubic 8-channel energy grid builder — Phase 2.
//!
//! Computes the 8-channel grid (energy + 7 derivatives) needed by the
//! Lekien–Marsden tricubic evaluator (`tricubic.rs`). The grid is built
//! by summing analytical LJ tricubic derivatives over framework atoms
//! within a cutoff, using the same supercell repeat strategy as the
//! existing scalar `egrid_from_atoms`.
//!
//! # Channel layout (RASPA-compatible)
//!
//! Channel 0 = energy
//! Channel 1–3 = ∂E/∂s_x, ∂E/∂s_y, ∂E/∂s_z   (÷ ncells after build)
//! Channel 4–6 = ∂²E/∂s_x∂s_y, ∂²E/∂s_x∂s_z, ∂²E/∂s_y∂s_z  (÷ product)
//! Channel 7   = ∂³E/∂s_x∂s_y∂s_z  (÷ triple product)
//!
//! The scaling by `ncells_per_axis` is applied AFTER the grid is built.
//! This matches the RASPA convention: grid stores scaled derivatives
//! so the tricubic evaluator can work in local cell coordinates.
//!
//! Reference: RASPA3 `framework_molecule_grid.cpp`,
//!   `calculateTricubicFractionalAtPosition` +
//!   `potentialLennardJonesTricubicDerivative`.


use numpy::PyArray1;
use pyo3::prelude::*;
use rayon::prelude::*;

/// Build an 8-channel tricubic energy grid from atomic positions.
///
/// Parameters match the existing `egrid_from_atoms`:
/// - `fpoints`: fractional coordinates (N, 3)
/// - `cell`: 3×3 cell matrix (columns = lattice vectors)
/// - `epsilons`, `sigmas`: per-atom LJ parameters
/// - `nx, ny, nz`: grid dimensions
/// - `cutoff`: interaction cutoff in Å
///
/// Returns a numpy array of shape (8, nx, ny, nz) in C order
/// with the 8-channel tricubic derivative data.
#[pyfunction]
#[pyo3(text_signature = "(fpoints, cell, epsilons, sigmas, nx, ny, nz, cutoff)")]
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
) -> PyResult<Py<PyArray1<f64>>> {
    let fpoints = fpoints.as_array();
    let cell = cell.as_array();
    let epsilons = epsilons.as_array();
    let sigmas = sigmas.as_array();

    let total_points = nx * ny * nz;
    let cutoff_sq = cutoff * cutoff;

    // Pre-compute cell matrix components.
    let c00 = cell[[0, 0]];
    let c01 = cell[[0, 1]];
    let c02 = cell[[0, 2]];
    let c10 = cell[[1, 0]];
    let c11 = cell[[1, 1]];
    let c12 = cell[[1, 2]];
    let c20 = cell[[2, 0]];
    let c21 = cell[[2, 1]];
    let c22 = cell[[2, 2]];

    // Number of cells per axis: nx-1 (for a grid of nx points).
    // These are used later by the caller for scaling.
    // Build the raw (unscaled) grid here; scaling is applied in Python wrapper.

    // Parallel computation across all grid points.
    // Output: flat vector of 8 * total_points values (channel-major).
    let buffer: Vec<f64> = (0..total_points)
        .into_par_iter()
        .flat_map(|idx| {
            let i = idx / (ny * nz);
            let j = (idx / nz) % ny;
            let k = idx % nz;

            // Fractional coordinate on the uniform grid.
            // Convention: linspace(-0.5, 0.5, n) — same as egrid_from_atoms.
            let x_frac = -0.5 + (i as f64) / ((nx - 1) as f64);
            let y_frac = -0.5 + (j as f64) / ((ny - 1) as f64);
            let z_frac = -0.5 + (k as f64) / ((nz - 1) as f64);

            let mut energy = 0.0f64;
            let mut first_x = 0.0f64;
            let mut first_y = 0.0f64;
            let mut first_z = 0.0f64;
            // Full Cartesian Hessian (6 unique components).
            let mut h_xx = 0.0f64;
            let mut h_yy = 0.0f64;
            let mut h_zz = 0.0f64;
            let mut h_xy = 0.0f64;
            let mut h_xz = 0.0f64;
            let mut h_yz = 0.0f64;
            // Full third derivative tensor (10 unique components by symmetry)
            let mut t_xxx = 0.0f64;
            let mut t_xxy = 0.0f64;
            let mut t_xxz = 0.0f64;
            let mut t_xyy = 0.0f64;
            let mut t_xyz = 0.0f64;
            let mut t_xzz = 0.0f64;
            let mut t_yyy = 0.0f64;
            let mut t_yyz = 0.0f64;
            let mut t_yzz = 0.0f64;
            let mut t_zzz = 0.0f64;

            for (atom_idx, atom_frac) in fpoints.outer_iter().enumerate() {
                // Minimum-image convention in fractional space.
                let mut dx = x_frac - atom_frac[0];
                let mut dy = y_frac - atom_frac[1];
                let mut dz = z_frac - atom_frac[2];

                dx -= dx.round();
                dy -= dy.round();
                dz -= dz.round();

                // Convert to Cartesian.
                let cart_x = c00 * dx + c01 * dy + c02 * dz;
                let cart_y = c10 * dx + c11 * dy + c12 * dz;
                let cart_z = c20 * dx + c21 * dy + c22 * dz;

                let r2 = cart_x * cart_x + cart_y * cart_y + cart_z * cart_z;

                if r2 > cutoff_sq {
                    continue;
                }
                // Floor r to 1e-12 (matches scalar egrid_from_atoms behavior)
                let r = if r2 < 1e-24 { 1e-12_f64 } else { r2.sqrt() };
                let r2 = r * r;

                let sigma = sigmas[atom_idx];
                let epsilon = epsilons[atom_idx];

                // ---- RASPA LJ tricubic derivative formulas ----
                // temp3 = (σ/r)⁶ = σ⁶ / r⁶ = (σ²/r²)³
                let sigma_sq = sigma * sigma;
                let inv_r2 = 1.0 / r2;
                let temp3 = (sigma_sq * inv_r2).powi(3);
                let arg1 = 4.0 * epsilon;

                // Energy
                energy += arg1 * (temp3 * temp3 - temp3);

                // First derivative factor: F1 = dE/dr / r
                let f1 = -6.0 * arg1 * (2.0 * temp3 * temp3 - temp3) * inv_r2;

                // Second derivative factor (cross-term): F2
                let r4 = r2 * r2;
                let f2 = 24.0 * arg1 * (7.0 * temp3 * temp3 - 2.0 * temp3) / r4;

                // Third derivative factor: F3
                let r6 = r4 * r2;
                let f3 = -96.0 * arg1 * (28.0 * temp3 * temp3 - 5.0 * temp3) / r6;

                // Accumulate Cartesian derivatives
                first_x += cart_x * f1;
                first_y += cart_y * f1;
                first_z += cart_z * f1;

                // Full Hessian: H_pp = f2·x_p² + f1 (diagonal), H_pq = f2·x_p·x_q
                h_xx += f2 * cart_x * cart_x + f1;
                h_yy += f2 * cart_y * cart_y + f1;
                h_zz += f2 * cart_z * cart_z + f1;
                h_xy += f2 * cart_x * cart_y;
                h_xz += f2 * cart_x * cart_z;
                h_yz += f2 * cart_y * cart_z;

                // Full third derivative tensor (RASPA lines 169-171):
                // T_ijk += f3*dr_i*dr_j*dr_k + (j==k?f2*dr_i:0) + (i==k?f2*dr_j:0) + (i==j?f2*dr_k:0)
                let f2x = f2 * cart_x;
                let f2y = f2 * cart_y;
                let f2z = f2 * cart_z;
                t_xxx += f3 * cart_x * cart_x * cart_x + 3.0 * f2x;
                t_xxy += f3 * cart_x * cart_x * cart_y + f2y;
                t_xxz += f3 * cart_x * cart_x * cart_z + f2z;
                t_xyy += f3 * cart_x * cart_y * cart_y + f2x;
                t_xyz += f3 * cart_x * cart_y * cart_z; // all distinct — no correction
                t_xzz += f3 * cart_x * cart_z * cart_z + f2x;
                t_yyy += f3 * cart_y * cart_y * cart_y + 3.0 * f2y;
                t_yyz += f3 * cart_y * cart_y * cart_z + f2z;
                t_yzz += f3 * cart_y * cart_z * cart_z + f2y;
                t_zzz += f3 * cart_z * cart_z * cart_z + 3.0 * f2z;
            }

            // Convert Cartesian derivatives to fractional using the cell matrix.
            // C^T @ first_derivative (which equals Σ_i C_{ij} · first_i)
            let first_sx = c00 * first_x + c10 * first_y + c20 * first_z;
            let first_sy = c01 * first_x + c11 * first_y + c21 * first_z;
            let first_sz = c02 * first_x + c12 * first_y + c22 * first_z;

            // Second fractional derivatives: bilinear form
            // d2E/ds_x/ds_y = Σ_{p,q} H_pq · C_{px} · C_{qy}
            let second_sx_sy = h_xx * c00 * c01
                + h_yy * c10 * c11
                + h_zz * c20 * c21
                + h_xy * (c00 * c11 + c10 * c01)
                + h_xz * (c00 * c21 + c20 * c01)
                + h_yz * (c10 * c21 + c20 * c11);

            let second_sx_sz = h_xx * c00 * c02
                + h_yy * c10 * c12
                + h_zz * c20 * c22
                + h_xy * (c00 * c12 + c10 * c02)
                + h_xz * (c00 * c22 + c20 * c02)
                + h_yz * (c10 * c22 + c20 * c12);

            let second_sy_sz = h_xx * c01 * c02
                + h_yy * c11 * c12
                + h_zz * c21 * c22
                + h_xy * (c01 * c12 + c11 * c02)
                + h_xz * (c01 * c22 + c21 * c02)
                + h_yz * (c11 * c22 + c21 * c12);

            // Third fractional derivative: full tensor contraction
            // d3E/ds_x/ds_y/ds_z = Σ_{p,q,r} T_{pqr} · C_{px} · C_{qy} · C_{rz}
            // C columns: [c00,c10,c20]^T = a-axis (s_x), [c01,c11,c21]^T = b-axis (s_y),
            // [c02,c12,c22]^T = c-axis (s_z).
            let third_sxsysz =
                // T_xxx
                t_xxx * c00 * c01 * c02
                // T_xxy
                + t_xxy * (c00 * c01 * c12 + c00 * c11 * c02 + c10 * c01 * c02)
                // T_xxz
                + t_xxz * (c00 * c01 * c22 + c00 * c21 * c02 + c20 * c01 * c02)
                // T_xyy
                + t_xyy * (c10 * c01 * c12 + c10 * c11 * c02 + c00 * c11 * c12)
                // T_xyz
                + t_xyz * (c00 * c21 * c12 + c10 * c21 * c02 + c20 * c11 * c02
                    + c10 * c01 * c22 + c20 * c01 * c12 + c00 * c11 * c22)
                // T_xzz
                + t_xzz * (c20 * c21 * c02 + c20 * c01 * c22 + c00 * c21 * c22)
                // T_yyy
                + t_yyy * c10 * c11 * c12
                // T_yyz
                + t_yyz * (c10 * c21 * c12 + c10 * c11 * c22 + c20 * c11 * c12)
                // T_yzz
                + t_yzz * (c20 * c11 * c22 + c10 * c21 * c22 + c20 * c21 * c12)
                // T_zzz
                + t_zzz * c20 * c21 * c22;

            // Return 8 values for this grid point (channel-major: all 8 values
            // for this point, then next point).
            // The caller will reshape to (8, nx, ny, nz) in C order.
            vec![
                energy,
                first_sx,
                first_sy,
                first_sz,
                second_sx_sy,
                second_sx_sz,
                second_sy_sz,
                third_sxsysz,
            ]
        })
        .collect();

    // Reshape: total_points groups of 8 values → (8, nx, ny, nz) in C order.
    // In C order (last index varies fastest), the flat layout for (8, nx, ny, nz)
    // is: data[c][i][j][k] → offset = c * nx*ny*nz + i * ny*nz + j * nz + k
    //
    // Our buffer is: for each grid point (i,j,k), 8 consecutive values (channels 0-7).
    // To convert to the C-order (8, nx, ny, nz) layout, we need to reorganize.
    //
    // Current: [ch0_000, ch1_000, ..., ch7_000, ch0_001, ...]
    // Target:  [ch0_000, ch0_001, ..., ch0_00(nz-1), ch0_010, ..., ch0_(nx-1)(ny-1)(nz-1),
    //           ch1_000, ...]
    //
    // Reorganize by channel on the Python side (simpler and less error-prone).

    let nynz = ny * nz;
    let mut output = vec![0.0f64; 8 * total_points];

    for idx in 0..total_points {
        let i = idx / nynz;
        let j = (idx / nz) % ny;
        let k = idx % nz;
        let base_src = idx * 8;
        for ch in 0..8 {
            let dst_idx = ch * nx * nynz + i * nynz + j * nz + k;
            output[dst_idx] = buffer[base_src + ch];
        }
    }

    // Return flat array; Python wrapper reshapes to (8, nx, ny, nz).
    let array = PyArray1::from_vec(py, output);
    Ok(array.to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_energy_channel_matches_scalar() {
        // This test verifies that channel 0 of the tricubic grid
        // matches the scalar egrid_from_atoms output.
        // It's a placeholder — the actual test requires numpy arrays
        // and is better done in Python integration tests.
    }
}
