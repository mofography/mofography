"""
Compare watershed backends (skimage, rust_bucket, rust_det) on MOF distance fields.

This script expects CIFs in ~/repos/mofography_testing and runs the
periodic segmentation pipeline with different watershed backends,
printing simple voxel disagreement statistics.
"""

from pathlib import Path

import numpy as np
import ase.io

import mofography as mgr
from mofography.rust_fun import dgrid_from_atoms_aabb_rust, do_periodic_segmentation_rust


TEST_ROOT = Path.home() / "repos" / "mofography_testing"


def _load_cif(name: str):
    path = TEST_ROOT / name
    if not path.is_file():
        print(f"[SKIP] CIF not found: {path}")
        return None
    return ase.io.read(path)


def _run_for_material(cif_name: str, spacing: float = 0.5):
    mat_atoms = _load_cif(cif_name)
    if mat_atoms is None:
        return

    print(f"\n=== {cif_name} ===")
    dgrid = dgrid_from_atoms_aabb_rust(mat_atoms, spacing=spacing)

    common_kwargs = dict(
        mask_thickness=0.5,
        h=0.2,
        spacing=spacing,
        threshold_abs=1.5,
        min_distance=1.5,
        compute_connections=False,
        return_dgrid=False,
    )

    res_sk = do_periodic_segmentation_rust(
        dgrid, mat_atoms, watershed_backend="skimage", **common_kwargs
    )
    res_bucket = do_periodic_segmentation_rust(
        dgrid, mat_atoms, watershed_backend="rust_bucket", **common_kwargs
    )
    res_det = do_periodic_segmentation_rust(
        dgrid, mat_atoms, watershed_backend="rust_det", **common_kwargs
    )

    reg_sk = res_sk["regions_pbc"].astype(np.int64)
    reg_bucket = res_bucket["regions_pbc"].astype(np.int64)
    reg_det = res_det["regions_pbc"].astype(np.int64)

    # Simple voxel-wise disagreement counts
    neq_bucket = np.count_nonzero(reg_sk != reg_bucket)
    neq_det = np.count_nonzero(reg_sk != reg_det)
    total = reg_sk.size

    print(f"voxels: {total}")
    print(f"skimage vs rust_bucket: {neq_bucket} differing voxels ({neq_bucket/total:.3%})")
    print(f"skimage vs rust_det   : {neq_det} differing voxels ({neq_det/total:.3%})")


def main():
    for cif in ["irmof-1.cif", "pim-0_4.cif"]:
        _run_for_material(cif)


if __name__ == "__main__":
    main()

