#!/usr/bin/env python
"""
Check consistency between pocket radii and window radii for a given CIF.

Usage
-----

    python devtools/check_pim_windows.py /path/to/pim-0_4.cif

This script:
- builds a distance grid using the Rust backend and TRAPPE-style radii,
- runs do_periodic_segmentation_rust with compute_connections=True,
- builds the directed pore graph, and
- reports any edges where window_radius > min(pocket_radius_i, pocket_radius_j).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read

import mofography as mgr
from mofography.winding_percolation import (
    make_directed_graph_from_segmentation_results,
)


def build_trappe_radii() -> pd.Series:
    """
    Return TRAPPE diameters as a pandas Series, then convert to radii on demand.

    This table is copied from the manuscript repro scripts so that the
    segmentation parameters for pim-0_4 are reproducible here.
    """

    from pandas import Series

    trappe_dia = Series(
        {
            "C1": 4.65,
            "C2": 3.95,
            "C3": 3.88,
            "C4": 3.695,
            "C5": 3.95,
            "C6": 4.65,
            "C7": 3.95,
            "O1": 2.599999,
            "C8": 3.7,
            "C9": 3.695002,
            "C10": 3.879998,
            "C11": 6.400003,
            "C12": 3.950004,
            "C13": 3.749998,
            "O2": 2.599999,
            "C14": 3.720001,
            "O3": 3.050001,
            "O4": 3.020004,
            "H1": 0.0,
            "C15": 3.7,
            "O5": 2.599999,
            "C16": 3.7,
            "C17": 3.695002,
            "C18": 3.879998,
            "C19": 6.400003,
            "C20": 3.950004,
            "C21": 3.749998,
            "C22": 3.7,
            "O6": 2.599999,
            "O8": 2.599999,
            "C23": 3.7,
            "C24": 3.695002,
            "C25": 3.879998,
            "C26": 6.400003,
            "C27": 3.950004,
            "C28": 3.749998,
            "C29": 3.7,
            "C30": 3.55,
            "N1": 2.949997,
            "O7": 2.599999,
            "C33": 3.879998,
            "C37": 3.7,
            "C31": 3.7,
            "O10": 2.599999,
            "C32": 3.695002,
            "C34": 6.400003,
            "C35": 3.950004,
            "C36": 3.749998,
            "F1": 3.080817,
            "O9": 2.599999,
            "O11": 2.599999,
            "C38": 3.7,
            "C39": 3.695002,
            "C40": 3.879998,
            "C41": 6.400003,
            "C42": 3.950004,
            "C43": 3.749998,
            "O12": 2.599999,
            "C44": 3.7,
            "O14": 2.599999,
            "C45": 3.7,
            "C46": 3.695002,
            "C47": 3.879998,
            "C48": 6.400003,
            "C49": 3.950004,
            "C491": 3.7,
            "C50": 3.55,
            "N2": 2.949997,
            "O13": 2.599999,
            "N3": 3.200001,
            "C51": 3.7,
            "C57": 3.7,
            "O16": 2.599999,
            "C52": 3.695002,
            "C53": 3.879998,
            "C54": 6.400003,
            "C55": 3.950004,
            "C56": 3.749998,
            "O15": 2.599999,
            "O18": 2.599999,
            "C58": 3.7,
            "C59": 3.695002,
            "C60": 3.879998,
            "C61": 6.400003,
            "C62": 3.950004,
            "C63": 3.749998,
            "C65": 3.7,
            "C66": 3.55,
            "N4": 2.949997,
            "O17": 2.599999,
            "O20": 2.599999,
            "C69": 3.879998,
            "C68": 3.695002,
            "C70": 6.400003,
            "C71": 3.950004,
            "C72": 3.749998,
            "C67": 3.7,
            "C75": 3.55,
            "N5": 2.949997,
            "C74": 3.7,
            "O19": 2.599999,
            "C78": 3.879998,
            "C76": 3.7,
            "O22": 2.599999,
            "O21": 2.599999,
            "C77": 3.695002,
            "C79": 6.400003,
            "C80": 3.950004,
            "C84": 3.55,
            "N6": 2.949997,
            "C81": 3.749998,
            "C83": 3.7,
            "C87": 3.879998,
            "C92": 3.7,
            "C85": 3.7,
            "O24": 2.599999,
            "C86": 3.695002,
            "C88": 6.400003,
            "C89": 3.950004,
            "C93": 3.55,
            "N7": 2.949997,
            "Si1": 3.804138,
            "C90": 3.749998,
            "O23": 2.599999,
            "C100": 3.7,
            "C101": 3.695002,
            "C102": 3.879998,
            "C103": 6.400003,
            "C104": 3.950004,
            "S2": 3.55,
            "C105": 3.749998,
            "C106": 3.55,
            "N8": 2.949997,
            "O26": 3.050001,
            "C94": 3.7,
            "C95": 3.695002,
            "C96": 3.879998,
            "C97": 6.400003,
            "C98": 3.950004,
            "S1": 3.55,
            "C99": 3.749998,
            "O25": 3.050001,
            "C113": 3.7,
            "C114": 3.695002,
            "C115": 3.879998,
            "C116": 6.400003,
            "C117": 3.950004,
            "S4": 3.55,
            "C118": 3.749998,
            "C119": 3.55,
            "N9": 2.949997,
            "C107": 3.7,
            "C108": 3.695002,
            "C109": 3.879998,
            "C110": 6.400003,
            "C111": 3.950004,
            "S3": 3.55,
            "C112": 3.749998,
        }
    )
    return trappe_dia


def run_single_segmentation(
    pim_path: Path,
    spacing: float = 0.25,
    mask_thickness: float = 0.0,
) -> tuple[dict, dict]:
    """
    Run periodic segmentation on a single CIF and return (results, metrics).
    """

    params = dict(
        spacing=spacing,
        h=0.5,
        threshold_abs=1.0,
        preserve_threshold=103.0,
        smooth_threshold=106.0,
        max_smooth_radius=2.0,
        extra_padding=2.0,
        min_distance=1.0,
        min_region_size=0.001,
        minimum_window_separation=1.5,
    )

    mat_atoms = read(str(pim_path), store_tags=True)
    site_labels = mat_atoms.info["_atom_site_label"]
    trappe_dia = build_trappe_radii()
    radii = trappe_dia[site_labels].values / 2.0

    dgrid = mgr.dgrid_from_atoms_rust(
        mat_atoms, spacing=params["spacing"], radii=radii
    )

    results = mgr.do_periodic_segmentation_rust(
        dgrid,
        mat_atoms,
        mask_thickness=mask_thickness,
        h=params["h"],
        spacing=params["spacing"],
        return_dgrid=True,
        threshold_abs=params["threshold_abs"],
        min_distance=params["min_distance"],
        return_padded_regions=True,
        preserve_threshold=params["preserve_threshold"],
        smooth_threshold=params["smooth_threshold"],
        max_smooth_radius=params["max_smooth_radius"],
        extra_padding=params["extra_padding"],
        compute_connections=True,
        min_region_size=params["min_region_size"],
        minimum_window_separation=params["minimum_window_separation"],
    )

    graph = make_directed_graph_from_segmentation_results(results, dgrid=results["dgrid"])

    metrics = dict(
        n_regions=int(results["regions"].max()),
        n_maxima=int(results["all_maxima"].shape[0]),
        n_connections=len(results["connections"]),
        plr=None,
        pld=None,
    )

    try:
        base_labels = mgr.build_base_region_labels_from_results(results)
        plr = mgr.find_pore_limiting_radius_winding(
            results, graph=graph, base_region_labels=base_labels, include_dangling=True
        )
        metrics["plr"] = plr
        metrics["pld"] = None if plr is None else 2.0 * plr
    except Exception:
        # PLR is not essential for the window vs pocket-radius check
        pass

    return results, graph, metrics


def find_radius_violations(
    graph,
    tol: float = 1e-3,
) -> list[dict]:
    """
    Return a list of edges where window_radius > min(pocket_radius_i, pocket_radius_j) + tol.
    """

    violations: list[dict] = []
    for u, v, data in graph.edges(data=True):
        w = float(data.get("window_radius", 0.0))
        ru = float(graph.nodes[u].get("pocket_radius", 0.0))
        rv = float(graph.nodes[v].get("pocket_radius", 0.0))
        limit = min(ru, rv) + tol
        if w > limit:
            violations.append(
                dict(
                    u=int(u),
                    v=int(v),
                    window_radius=w,
                    pocket_radius_u=ru,
                    pocket_radius_v=rv,
                    limit=limit,
                )
            )
    return violations


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("Usage: python devtools/check_pim_windows.py /path/to/pim-0_4.cif [mask_thickness]")
        return 1

    pim_path = Path(argv[1]).expanduser()
    if not pim_path.exists():
        print(f"Error: CIF not found at {pim_path}")
        return 1

    mask_thickness = float(argv[2]) if len(argv) > 2 else 0.0

    print(f"[info] CIF: {pim_path}")
    print(f"[info] mask_thickness = {mask_thickness}")

    results, graph, metrics = run_single_segmentation(
        pim_path, spacing=0.25, mask_thickness=mask_thickness
    )

    print(
        f"[info] regions={metrics['n_regions']}, "
        f"maxima={metrics['n_maxima']}, "
        f"connections={metrics['n_connections']}, "
        f"plr={metrics['plr']}, pld={metrics['pld']}"
    )

    violations = find_radius_violations(graph, tol=1e-3)
    if not violations:
        print("[ok] No window-radius > pocket-radius violations found.")
        return 0

    print(
        f"[warn] Found {len(violations)} edges where "
        "window_radius > min(pocket_radius_i, pocket_radius_j) + tol"
    )
    # Sort for stable output and highlight node 20 if present
    violations_sorted = sorted(
        violations, key=lambda d: (min(d["u"], d["v"]), max(d["u"], d["v"]), -d["window_radius"])
    )
    for v in violations_sorted:
        u = v["u"]
        w = v["v"]
        mark = "  "
        if 20 in (u, w):
            mark = "**"
        print(
            f"{mark}edge {u}-{w}: "
            f"window_radius={v['window_radius']:.3f} Å, "
            f"pocket_u={v['pocket_radius_u']:.3f} Å, "
            f"pocket_v={v['pocket_radius_v']:.3f} Å"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

