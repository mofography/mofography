import hashlib
import os
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read

import mofography as mgr
from mofography import rust_fun as rust_mod
from mofography import winding_percolation as wp
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def _hash_array(arr: np.ndarray) -> str:
    """Return a short hash for an integer array (order-insensitive for 1D labels)."""
    arr = np.asarray(arr)
    if arr.size == 0:
        data = np.empty((0,), dtype=np.int64)
    elif arr.ndim == 1:
        data = np.sort(arr.astype(np.int64))
    else:
        flat = arr.reshape(arr.shape[0], -1)
        order = np.lexsort(flat.T)
        data = flat[order].astype(np.int64)
    return hashlib.md5(data.tobytes()).hexdigest()[:8]


def build_trappe_radii():
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


BACKEND = os.getenv("WATERSHED_BACKEND", "rust_det_parallel")


def run_single_segmentation(
    dgrid,
    mat_atoms,
    mask_thickness: float,
    params: dict,
) -> tuple[dict, np.ndarray]:
    results = rust_mod.do_periodic_segmentation_rust(
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
        watershed_backend=BACKEND,
    )

    regions = results["regions"]
    all_maxima = results["all_maxima"]
    maxima_pbc = results["maxima_pbc"]
    connections = results["connections"]

    graph = wp.make_directed_graph_from_segmentation_results(results)
    base_labels = wp.build_base_region_labels_from_results(results)
    plr = wp.find_pore_limiting_radius_winding(
        results, graph=graph, base_region_labels=base_labels, include_dangling=True
    )
    pld = None if plr is None else 2 * plr

    if plr is not None:
        axis_clusters, _, _ = wp.percolation_all_axes_winding(
            results,
            threshold=plr,
            graph=graph,
            base_region_labels=base_labels,
            include_dangling=True,
        )
        percolation_regions = np.unique(
            np.array(
                axis_clusters["x"] + axis_clusters["y"] + axis_clusters["z"]
            ).flatten()
        )
    else:
        percolation_regions = np.array([], dtype=int)

    # Basic metrics and hashes for reproducibility checks
    maxima_hash = _hash_array(all_maxima.reshape(all_maxima.shape[0], -1))
    maxima_pbc_hash = _hash_array(maxima_pbc.reshape(maxima_pbc.shape[0], -1))
    regions_hash = _hash_array(np.unique(regions))
    conn_pairs = np.array([c[0] for c in connections], dtype=int) if connections else np.zeros((0, 2), dtype=int)
    connections_hash = _hash_array(conn_pairs)
    perc_hash = _hash_array(percolation_regions)

    metrics = dict(
        mask_thickness=mask_thickness,
        n_regions=int(regions.max()),
        n_maxima=int(all_maxima.shape[0]),
        n_maxima_pbc=int(maxima_pbc.shape[0]),
        n_connections=len(connections),
        plr=plr,
        pld=pld,
        n_percolation_regions=int(percolation_regions.size),
        maxima_hash=maxima_hash,
        maxima_pbc_hash=maxima_pbc_hash,
        regions_hash=regions_hash,
        connections_hash=connections_hash,
        percolation_hash=perc_hash,
    )
    return metrics, regions


def main():
    base = Path(__file__).resolve().parent
    pim_path = (base / "pim-0_4.cif").expanduser()
    if not pim_path.exists():
        raise FileNotFoundError(f"Cannot find CIF at {pim_path}")

    # Parameters matching the existing debug org for pim-0_4
    params = dict(
        spacing=0.5,
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

    mat_atoms = read(pim_path, store_tags=True)
    site_labels = mat_atoms.info["_atom_site_label"]
    trappe_dia = build_trappe_radii()
    radii = trappe_dia[site_labels].values / 2.0

    dgrid = rust_mod.dgrid_from_atoms_rust(mat_atoms, spacing=params["spacing"], radii=radii)

    mask_values = [0.0, 0.25]
    n_trials = 3

    records = []
    for mask_thickness in mask_values:
        regions_ref = None
        boundary_ref = None
        for trial in range(1, n_trials + 1):
            metrics, regions = run_single_segmentation(
                dgrid, mat_atoms, mask_thickness, params
            )
            # quantify voxel-level variability vs first trial at this mask_thickness
            if regions_ref is None:
                regions_ref = regions.copy()
                # define boundary voxels on the reference labeling
                boundary_ref = np.zeros_like(regions_ref, dtype=bool)
                for axis in range(3):
                    sl1 = [slice(None)] * 3
                    sl2 = [slice(None)] * 3
                    sl1[axis] = slice(0, -1)
                    sl2[axis] = slice(1, None)
                    neq = regions_ref[tuple(sl1)] != regions_ref[tuple(sl2)]
                    boundary_ref[tuple(sl1)] |= neq
                    boundary_ref[tuple(sl2)] |= neq
                frac_changed = 0.0
                frac_changed_boundary = 0.0
            else:
                # compare only over voxels that are labeled in at least one run
                union_mask = (regions_ref != 0) | (regions != 0)
                diff = regions != regions_ref
                if union_mask.any():
                    frac_changed = float(diff[union_mask].mean())
                    # boundary voxels subset
                    boundary_union = union_mask & boundary_ref
                    if boundary_union.any():
                        frac_changed_boundary = float(
                            diff[boundary_union].mean()
                        )
                    else:
                        frac_changed_boundary = 0.0
                else:
                    frac_changed = 0.0
                    frac_changed_boundary = 0.0
            metrics["frac_voxels_label_changed"] = frac_changed
            metrics["frac_boundary_voxels_label_changed"] = frac_changed_boundary
            metrics["trial"] = trial
            records.append(metrics)

    df = pd.DataFrame.from_records(records)
    df = df[
        [
            "mask_thickness",
            "trial",
            "frac_voxels_label_changed",
            "frac_boundary_voxels_label_changed",
            "n_regions",
            "n_maxima",
            "n_maxima_pbc",
            "n_connections",
            "plr",
            "pld",
            "n_percolation_regions",
            "maxima_hash",
            "maxima_pbc_hash",
            "regions_hash",
            "connections_hash",
            "percolation_hash",
        ]
    ]

    out_csv = base / "repro_segmentation_results.csv"
    df.to_csv(out_csv, index=False)

    # Also emit a simple org-mode table for quick inspection
    out_org = base / "repro_segmentation_results.org"
    with out_org.open("w") as f:
        f.write("* Segmentation reproducibility for pim-0_4\n")
        f.write("  Generated by repro_segmentation_tests.py\n\n")
        f.write("| mask_thickness | trial | frac_voxels_label_changed | frac_boundary_voxels_label_changed | n_regions | n_maxima | n_maxima_pbc | n_connections | plr | pld | n_percolation_regions | maxima_hash | maxima_pbc_hash | regions_hash | connections_hash | percolation_hash |\n")
        f.write("|-\n")
        for _, row in df.iterrows():
            f.write(
                f"| {row['mask_thickness']:.2f} | {int(row['trial'])} | "
                f"{row['frac_voxels_label_changed']:.6f} | {row['frac_boundary_voxels_label_changed']:.6f} | "
                f"{int(row['n_regions'])} | {int(row['n_maxima'])} | {int(row['n_maxima_pbc'])} | "
                f"{int(row['n_connections'])} | {row['plr'] if pd.notna(row['plr']) else 'NA'} | "
                f"{row['pld'] if pd.notna(row['pld']) else 'NA'} | "
                f"{int(row['n_percolation_regions'])} | "
                f"{row['maxima_hash']} | {row['maxima_pbc_hash']} | {row['regions_hash']} | "
                f"{row['connections_hash']} | {row['percolation_hash']} |\n"
            )

    print(f"Wrote {out_csv}")
    print(f"Wrote {out_org}")


if __name__ == "__main__":
    main()

