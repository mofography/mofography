//! Shared MIC, Cartesian derivative accumulation, and fractional contraction
//! for 8-channel tricubic grid builders (LJ, Coulomb).

/// Precomputed cell-matrix products for Cartesian → fractional contraction.
pub struct CellContraction {
    pub c00: f64,
    pub c01: f64,
    pub c02: f64,
    pub c10: f64,
    pub c11: f64,
    pub c12: f64,
    pub c20: f64,
    pub c21: f64,
    pub c22: f64,
    c2_xx_xy: f64,
    c2_xx_xz: f64,
    c2_yy_xy: f64,
    c2_zz_xz: f64,
    c2_xy_cross: f64,
    c2_xz_cross: f64,
    c2_yz_cross: f64,
    cp_xy_xz: f64,
    cp_xy_yz: f64,
    cp_xx_yz: f64,
    cp_xx_zz: f64,
    cp_yy_xz: f64,
}

impl CellContraction {
    pub fn from_cell_matrix(
        c00: f64,
        c01: f64,
        c02: f64,
        c10: f64,
        c11: f64,
        c12: f64,
        c20: f64,
        c21: f64,
        c22: f64,
    ) -> Self {
        Self {
            c00,
            c01,
            c02,
            c10,
            c11,
            c12,
            c20,
            c21,
            c22,
            c2_xx_xy: c00 * c01,
            c2_xx_xz: c00 * c02,
            c2_yy_xy: c10 * c11,
            c2_zz_xz: c20 * c22,
            c2_xy_cross: c00 * c11 + c10 * c01,
            c2_xz_cross: c00 * c21 + c20 * c01,
            c2_yz_cross: c10 * c21 + c20 * c11,
            cp_xy_xz: c00 * c21 + c20 * c01,
            cp_xy_yz: c10 * c21 + c20 * c11,
            cp_xx_yz: c00 * c12 + c10 * c02,
            cp_xx_zz: c00 * c22 + c20 * c02,
            cp_yy_xz: c01 * c22 + c21 * c02,
        }
    }
}

/// Cartesian derivative tensors accumulated from radial f1, f2, f3.
#[derive(Clone, Copy, Default)]
pub struct CartesianDerivatives {
    pub energy: f64,
    pub first_x: f64,
    pub first_y: f64,
    pub first_z: f64,
    pub h_xx: f64,
    pub h_yy: f64,
    pub h_zz: f64,
    pub h_xy: f64,
    pub h_xz: f64,
    pub h_yz: f64,
    pub t_xxx: f64,
    pub t_xxy: f64,
    pub t_xxz: f64,
    pub t_xyy: f64,
    pub t_xyz: f64,
    pub t_xzz: f64,
    pub t_yyy: f64,
    pub t_yyz: f64,
    pub t_yzz: f64,
    pub t_zzz: f64,
}

impl CartesianDerivatives {
    pub fn accumulate_radial(&mut self, dx: f64, dy: f64, dz: f64, f1: f64, f2: f64, f3: f64) {
        self.first_x += dx * f1;
        self.first_y += dy * f1;
        self.first_z += dz * f1;

        self.h_xx += f2 * dx * dx + f1;
        self.h_yy += f2 * dy * dy + f1;
        self.h_zz += f2 * dz * dz + f1;
        self.h_xy += f2 * dx * dy;
        self.h_xz += f2 * dx * dz;
        self.h_yz += f2 * dy * dz;

        let f2x = f2 * dx;
        let f2y = f2 * dy;
        let f2z = f2 * dz;
        self.t_xxx += f3 * dx * dx * dx + 3.0 * f2x;
        self.t_xxy += f3 * dx * dx * dy + f2y;
        self.t_xxz += f3 * dx * dx * dz + f2z;
        self.t_xyy += f3 * dx * dy * dy + f2x;
        self.t_xyz += f3 * dx * dy * dz;
        self.t_xzz += f3 * dx * dz * dz + f2x;
        self.t_yyy += f3 * dy * dy * dy + 3.0 * f2y;
        self.t_yyz += f3 * dy * dy * dz + f2z;
        self.t_yzz += f3 * dy * dz * dz + f2y;
        self.t_zzz += f3 * dz * dz * dz + 3.0 * f2z;
    }

    pub fn to_fractional_channels(
        &self,
        cell: &CellContraction,
        nx: usize,
        ny: usize,
        nz: usize,
        extent_x: f64,
        extent_y: f64,
        extent_z: f64,
    ) -> [f64; 8] {
        let first_sx = cell.c00 * self.first_x + cell.c10 * self.first_y + cell.c20 * self.first_z;
        let first_sy = cell.c01 * self.first_x + cell.c11 * self.first_y + cell.c21 * self.first_z;
        let first_sz = cell.c02 * self.first_x + cell.c12 * self.first_y + cell.c22 * self.first_z;

        let second_sx_sy = self.h_xx * cell.c2_xx_xy
            + self.h_yy * cell.c2_yy_xy
            + self.h_zz * cell.c20 * cell.c21
            + self.h_xy * cell.c2_xy_cross
            + self.h_xz * cell.cp_xy_xz
            + self.h_yz * cell.cp_xy_yz;

        let second_sx_sz = self.h_xx * cell.c2_xx_xz
            + self.h_yy * cell.c10 * cell.c12
            + self.h_zz * cell.c2_zz_xz
            + self.h_xy * cell.cp_xx_yz
            + self.h_xz * cell.c2_xz_cross
            + self.h_yz * cell.cp_xx_zz;

        let second_sy_sz = self.h_xx * cell.c01 * cell.c02
            + self.h_yy * cell.c11 * cell.c12
            + self.h_zz * cell.c21 * cell.c22
            + self.h_xy * (cell.c01 * cell.c12 + cell.c11 * cell.c02)
            + self.h_xz * cell.cp_yy_xz
            + self.h_yz * cell.c2_yz_cross;

        let third_sxsysz = self.t_xxx * cell.c00 * cell.c01 * cell.c02
            + self.t_xxy
                * (cell.c00 * cell.c01 * cell.c12
                    + cell.c00 * cell.c11 * cell.c02
                    + cell.c10 * cell.c01 * cell.c02)
            + self.t_xxz
                * (cell.c00 * cell.c01 * cell.c22
                    + cell.c00 * cell.c21 * cell.c02
                    + cell.c20 * cell.c01 * cell.c02)
            + self.t_xyy
                * (cell.c10 * cell.c01 * cell.c12
                    + cell.c10 * cell.c11 * cell.c02
                    + cell.c00 * cell.c11 * cell.c12)
            + self.t_xyz
                * (cell.c00 * cell.c21 * cell.c12
                    + cell.c10 * cell.c21 * cell.c02
                    + cell.c20 * cell.c11 * cell.c02
                    + cell.c10 * cell.c01 * cell.c22
                    + cell.c20 * cell.c01 * cell.c12
                    + cell.c00 * cell.c11 * cell.c22)
            + self.t_xzz
                * (cell.c20 * cell.c21 * cell.c02
                    + cell.c20 * cell.c01 * cell.c22
                    + cell.c00 * cell.c21 * cell.c22)
            + self.t_yyy * cell.c10 * cell.c11 * cell.c12
            + self.t_yyz
                * (cell.c10 * cell.c21 * cell.c12
                    + cell.c10 * cell.c11 * cell.c22
                    + cell.c20 * cell.c11 * cell.c12)
            + self.t_yzz
                * (cell.c20 * cell.c11 * cell.c22
                    + cell.c10 * cell.c21 * cell.c22
                    + cell.c20 * cell.c21 * cell.c12)
            + self.t_zzz * cell.c20 * cell.c21 * cell.c22;

        let sx_scl = extent_x / ((nx - 1) as f64);
        let sy_scl = extent_y / ((ny - 1) as f64);
        let sz_scl = extent_z / ((nz - 1) as f64);

        [
            self.energy,
            first_sx * sx_scl,
            first_sy * sy_scl,
            first_sz * sz_scl,
            second_sx_sy * sx_scl * sy_scl,
            second_sx_sz * sx_scl * sz_scl,
            second_sy_sz * sy_scl * sz_scl,
            third_sxsysz * sx_scl * sy_scl * sz_scl,
        ]
    }
}

#[inline]
pub fn grid_fractional_coords(
    i: usize,
    j: usize,
    k: usize,
    nx: usize,
    ny: usize,
    nz: usize,
    extent_x: f64,
    extent_y: f64,
    extent_z: f64,
) -> (f64, f64, f64) {
    let denom_x = (nx - 1) as f64;
    let denom_y = (ny - 1) as f64;
    let denom_z = (nz - 1) as f64;
    (
        -0.5 + (i as f64) / denom_x * extent_x,
        -0.5 + (j as f64) / denom_y * extent_y,
        -0.5 + (k as f64) / denom_z * extent_z,
    )
}

#[inline]
pub fn min_image_frac_delta(
    x_frac: f64,
    y_frac: f64,
    z_frac: f64,
    atom_x: f64,
    atom_y: f64,
    atom_z: f64,
) -> (f64, f64, f64) {
    let mut dfx = x_frac - atom_x;
    let mut dfy = y_frac - atom_y;
    let mut dfz = z_frac - atom_z;
    dfx -= dfx.round();
    dfy -= dfy.round();
    dfz -= dfz.round();
    (dfx, dfy, dfz)
}

#[inline]
pub fn frac_delta_to_cart(
    cell: &CellContraction,
    dfx: f64,
    dfy: f64,
    dfz: f64,
) -> (f64, f64, f64) {
    let dx = cell.c00 * dfx + cell.c01 * dfy + cell.c02 * dfz;
    let dy = cell.c10 * dfx + cell.c11 * dfy + cell.c12 * dfz;
    let dz = cell.c20 * dfx + cell.c21 * dfy + cell.c22 * dfz;
    (dx, dy, dz)
}

pub fn reorganize_to_channel_major(buffer: Vec<f64>, nx: usize, ny: usize, nz: usize) -> Vec<f64> {
    let total_points = nx * ny * nz;
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
    output
}
