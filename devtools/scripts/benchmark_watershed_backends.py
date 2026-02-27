"""
Quick benchmarks for watershed backends:

- mgr_rust.seeded_watershed_det     (heap-based deterministic)
- mgr_rust.seeded_watershed_bucket  (parallel bucket)
- skimage.segmentation.watershed    (Python heap)

Runs on:
- Synthetic small 3D grids
- IRMOF-1 (P1-IRMOF-1.cif) at coarse spacing
- pim-0_4.cif at coarse spacing
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
from ase.io import read
from skimage.segmentation import watershed as sk_watershed

import mofography as mgr
from mofography import rust_fun as rust_mod


def _make_synthetic(shape=(64, 64, 64)):
    z, y, x = np.indices(shape)
    center1 = np.array(shape) // 3
    center2 = np.array(shape) * 2 // 3
    dist1 = np.sqrt(((z - center1[0]) ** 2 + (y - center1[1]) ** 2 + (x - center1[2]) ** 2))
    dist2 = np.sqrt(((z - center2[0]) ** 2 + (y - center2[1]) ** 2 + (x - center2[2]) ** 2))
    dgrid = np.minimum(dist1, dist2).astype(np.float32)
    seeds = np.zeros_like(dgrid, dtype=np.uint32)
    seeds[tuple(center1)] = 1
    seeds[tuple(center2)] = 2
    # For benchmarking, keep all voxels in the mask so seeds are always valid.
    mask = np.ones_like(dgrid, dtype=np.uint8)
    return dgrid, seeds, mask


def _time_call(label, fn, *args, repeat=3, **kwargs):
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn(*args, **kwargs)
        t1 = time.perf_counter()
        times.append(t1 - t0)
    times = np.array(times, dtype=float)
    print(f"{label:20s}  mean={times.mean():7.3f}s  min={times.min():7.3f}s  max={times.max():7.3f}s")


def bench_synthetic():
    print("=== Synthetic volume benchmarks ===")
    dgrid, seeds, mask = _make_synthetic()

    def sk_fn():
        sk_watershed(
            -dgrid,
            seeds.astype(np.int32),
            mask=mask.astype(bool),
            connectivity=2,
            compactness=0,
        )

    def det_fn():
        rust_mod.seeded_watershed_det_rust(dgrid, seeds, mask=mask)

    def bucket_fn():
        import mgr_rust

        mgr_rust.seeded_watershed_bucket(
            dgrid.astype(np.float32),
            seeds.astype(np.uint32),
            mask=mask.astype(np.uint8),
        )

    _time_call("skimage", sk_fn)
    _time_call("rust_det (heap)", det_fn)
    _time_call("rust_bucket", bucket_fn)


def bench_mof(cif_path: Path, name: str, spacing: float = 0.25):
    print(f"\n=== {name} (spacing={spacing}) benchmarks ===")
    mat_atoms = read(str(cif_path), store_tags=True)
    dgrid = rust_mod.dgrid_from_atoms_aabb_rust(mat_atoms, spacing=spacing)

    # For skimage, re-use markers from the Rust maxima detector to be fair
    res = rust_mod.do_periodic_segmentation_rust(
        dgrid,
        mat_atoms,
        mask_thickness=0.5,
        h=0.5,
        spacing=spacing,
        threshold_abs=1.0,
        min_distance=1.0,
        preserve_threshold=103.0,
        smooth_threshold=106.0,
        max_smooth_radius=2.0,
        extra_padding=2.0,
        min_region_size=0.001,
        minimum_window_separation=1.5,
        compute_connections=False,
        return_dgrid=False,
        watershed_backend="rust_bucket",
    )
    merged_maxima = res["maxima_pbc"]
    dt = dgrid.copy()
    # build markers in central (unpadded) grid, same as segmentation pipeline
    markers = np.zeros_like(dt, dtype=np.uint32)
    if merged_maxima.size > 0:
        coords = tuple(merged_maxima.T)
        markers[coords] = np.arange(1, len(merged_maxima) + 1, dtype=np.uint32)
    void_mask = (dt > 0.5).astype(np.uint8)

    def sk_fn():
        sk_watershed(
            -dt,
            markers.astype(np.int32),
            mask=void_mask.astype(bool),
            connectivity=2,
            compactness=0,
        )

    def det_fn():
        rust_mod.seeded_watershed_det_rust(dt, markers, mask=void_mask)

    def bucket_fn():
        import mgr_rust

        mgr_rust.seeded_watershed_bucket(
            dt.astype(np.float32),
            markers.astype(np.uint32),
            mask=void_mask.astype(np.uint8),
        )

    _time_call("skimage", sk_fn, repeat=3)
    _time_call("rust_det (heap)", det_fn, repeat=3)
    _time_call("rust_bucket", bucket_fn, repeat=3)


def main():
    bench_synthetic()

    repo_root = Path(__file__).resolve().parents[2]
    irmof = repo_root / "examples" / "TestCIFs" / "P1-IRMOF-1.cif"
    pim = Path.home() / "repos" / "mofography_testing" / "pim-0_4.cif"

    if irmof.is_file():
        bench_mof(irmof, "P1-IRMOF-1", spacing=0.25)
    else:
        print("Skipping IRMOF-1: CIF not found at", irmof)

    if pim.is_file():
        bench_mof(pim, "pim-0_4", spacing=0.25)
    else:
        print("Skipping pim-0_4: CIF not found at", pim)


if __name__ == "__main__":
    main()

