#!/usr/bin/env python3
"""Validate rust_det_metric against test_metric_watershed prototype on VIGNOQ_clean."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from ase.io import read

import mofography as mgr
import mgr_rust
import test_metric_watershed as tmw
from mofography.compute import calculate_padding_thickness, find_pore_limiting_radius
from mofography.rust_fun import (
    _compute_metric_inv_weights,
    _cell_matrix_for_metric_weights,
    dgrid_from_atoms_rust,
    do_periodic_segmentation_rust,
)

ROOT = Path(__file__).resolve().parent
VIGNOQ_CIF = ROOT / "examples" / "TestCIFs" / "VIGNOQ_clean.cif"
DEFAULT_HTML = ROOT / "vignoq_rust_det_metric.html"

SPACING = 0.5
MASK_THICKNESS = 0.05
H_PARAM = 0.2
EXTRA_PADDING_ANGSTROM = 1.0
PLD_TARGET = 1.378
PLD_TOL = 0.05
SELF_LOOP_MAX_TARGET = 1.38
SELF_LOOP_MAX_TOL = 0.15


def _proto_regions(dgrid, atoms, maxima_info, padding: int) -> np.ndarray:
    padded = np.pad(dgrid, padding, mode="wrap")
    void = np.pad(dgrid > MASK_THICKNESS, padding, mode="wrap")
    inv = _compute_metric_inv_weights(
        _cell_matrix_for_metric_weights(atoms), dgrid.shape, SPACING
    )
    markers = tmw.build_unmerged_markers_from_maxima_info(
        np.asarray(maxima_info, dtype=int), dgrid.shape, padding
    )
    labels = mgr_rust.seeded_watershed_det_metric(
        padded.astype(np.float32),
        markers.astype(np.uint32),
        inv,
        void.astype(np.uint8),
    )
    return tmw.crop_padded_labels(labels, padding, dgrid.shape)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cif", type=Path, default=VIGNOQ_CIF)
    parser.add_argument(
        "--html",
        type=Path,
        nargs="?",
        const=DEFAULT_HTML,
        default=None,
        help="Write Plotly HTML (default path if flag alone)",
    )
    args = parser.parse_args()

    atoms = read(str(args.cif.expanduser().resolve()))
    radii = mgr.default_vdw_radii_angstrom(atoms)
    dgrid = dgrid_from_atoms_rust(atoms, spacing=SPACING, radii=radii).astype(np.float32)
    padding = int(
        calculate_padding_thickness(dgrid, atoms, SPACING)
        + np.ceil(EXTRA_PADDING_ANGSTROM / SPACING)
    )

    plr = float(find_pore_limiting_radius(dgrid))
    print(f"PLD (dgrid percolation): {plr:.4f} A")

    bucket = do_periodic_segmentation_rust(
        dgrid,
        atoms,
        spacing=SPACING,
        mask_thickness=MASK_THICKNESS,
        h=H_PARAM,
        extra_padding=EXTRA_PADDING_ANGSTROM,
        compute_connections=False,
        watershed_backend="rust_bucket",
    )
    proto = _proto_regions(dgrid, atoms, bucket["maxima_info"], padding)

    metric = do_periodic_segmentation_rust(
        dgrid,
        atoms,
        spacing=SPACING,
        mask_thickness=MASK_THICKNESS,
        h=H_PARAM,
        extra_padding=EXTRA_PADDING_ANGSTROM,
        compute_connections=True,
        watershed_backend="rust_det_metric",
        _metric_bucket_scaffold=bucket,
    )
    regions = np.asarray(metric["regions"], dtype=np.uint32)

    if not np.array_equal(proto, regions):
        diff = int(np.sum(proto != regions))
        raise SystemExit(
            f"FAIL: rust_det_metric regions differ from Python prototype ({diff} voxels)"
        )
    print("PASS: rust_det_metric regions match Python prototype")

    # boundary stats on prototype labels
    padded = np.pad(dgrid, padding, mode="wrap")
    void = np.pad(dgrid > MASK_THICKNESS, padding, mode="wrap")
    inv = _compute_metric_inv_weights(
        _cell_matrix_for_metric_weights(atoms), dgrid.shape, SPACING
    )
    markers = tmw.build_unmerged_markers_from_maxima_info(
        np.asarray(bucket["maxima_info"], dtype=int), dgrid.shape, padding
    )
    labels = mgr_rust.seeded_watershed_det_metric(
        padded.astype(np.float32),
        markers.astype(np.uint32),
        inv,
        void.astype(np.uint8),
    )
    offsets_half = tmw.positive_half_offsets(tmw.neighbor_offsets_26())
    pbc_groups = [set(int(x) for x in g) for g in bucket["pbc_groups"]]
    pocket1, pocket2 = pbc_groups[0], pbc_groups[1]
    self_stats = tmw.summarize(
        tmw.boundary_values_for_label_sets(
            labels, padded, void, offsets_half, pocket1, None
        )
    )
    inter_stats = tmw.summarize(
        tmw.boundary_values_for_label_sets(
            labels, padded, void, offsets_half, pocket1, pocket2
        )
    )
    print(
        "Boundary self-loop: "
        f"max={self_stats['max']:.4f} A (target ~{SELF_LOOP_MAX_TARGET} A)"
    )
    print(
        "Boundary inter-pocket: "
        f"max={inter_stats['max']:.4f} A"
    )

    if abs(plr - PLD_TARGET) > PLD_TOL:
        raise SystemExit(
            f"FAIL: PLD {plr:.4f} A outside [{PLD_TARGET - PLD_TOL:.3f}, "
            f"{PLD_TARGET + PLD_TOL:.3f}] A"
        )

    conns = metric.get("connections") or []
    if conns:
        self_mins = []
        all_mins = []
        for c in conns:
            d = np.asarray(c[-1], dtype=float)
            if d.size == 0:
                continue
            b = float(np.min(d))
            all_mins.append(b)
            pair = [int(x) for x in c[0]]
            if pair[0] in pocket1 and pair[1] in pocket1:
                self_mins.append(b)
        print(
            "Connection bottlenecks: "
            f"max={max(all_mins):.4f} A, "
            f"self-loop max={max(self_mins):.4f} A (n={len(self_mins)})"
        )
        if self_mins and abs(max(self_mins) - SELF_LOOP_MAX_TARGET) > SELF_LOOP_MAX_TOL:
            raise SystemExit(
                f"FAIL: self-loop connection bottleneck {max(self_mins):.4f} A "
                f"outside ~{SELF_LOOP_MAX_TARGET} +/- {SELF_LOOP_MAX_TOL} A"
            )

    if abs(self_stats["max"] - SELF_LOOP_MAX_TARGET) > SELF_LOOP_MAX_TOL:
        print(
            "WARN: self-loop boundary max outside expected band "
            f"[{SELF_LOOP_MAX_TARGET - SELF_LOOP_MAX_TOL:.2f}, "
            f"{SELF_LOOP_MAX_TARGET + SELF_LOOP_MAX_TOL:.2f}] A"
        )

    if args.html is not None:
        tmw.write_visualization_html(
            out_path=args.html,
            mat_atoms=atoms,
            dgrid=dgrid,
            seg=metric,
            regions_standard=regions,
            regions_metric=regions,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
