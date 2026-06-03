"""
Build a Three.js HTML scene with framework, cell, regions, regions PBC,
maxima, pore graph, AND two transparent padded surface traces.

Usage (command-line):
    source ~/miniforge3/bin/activate mgr_env
    cd ~/repos/mofography
    python make_padded_threejs_scene.py --cif examples/TestCIFs/P1-IRMOF-1.cif
    python make_padded_threejs_scene.py --cif ~/repos/mgr_cli_py/cifs/P1-Cu-BTC.cif

Usage (Python API / workflow):
    from mofography.workflow import build_padded_threejs_scene
    html_path = build_padded_threejs_scene("path/to/file.cif", output_path="out.html")
"""
import argparse
import sys
from pathlib import Path
import numpy as np
from ase.io import read

# Ensure mofography is importable
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mofography as mgr
from mofography.movies import (
    MovieFigure,
    _make_trace,
    make_surface_trace,
)
from skimage.measure import marching_cubes
from scipy.interpolate import RegularGridInterpolator as RGI

try:
    import trimesh
    _HAS_TRIMESH = True
except Exception:
    _HAS_TRIMESH = False
    trimesh = None


# ---------------------------------------------------------------------------
# Helpers (same as previous version)
# ---------------------------------------------------------------------------

def _make_padded_surface_trace(
    dgrid,
    regions_padded,
    atoms,
    padding_thickness,
    isovalue=1.0,
    step_size=2,
    opacity=0.2,
    trace_id="surface_padded",
    name="Padded regions",
    cscale=None,
    select_labels=None,
):
    """Build a Three.js mesh trace for a padded segmentation surface."""
    pad = int(padding_thickness)

    padded_dgrid = np.pad(dgrid, pad, mode="wrap")
    assert padded_dgrid.shape == regions_padded.shape

    verts, faces, _normals, _vals = marching_cubes(
        padded_dgrid, level=isovalue, step_size=max(step_size, 1)
    )
    if len(faces) == 0:
        return {
            "positions": [],
            "indices": [],
            "colors": [],
            "n_vertices": 0,
            "n_triangles": 0,
        }

    if _HAS_TRIMESH:
        mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
        try:
            trimesh.smoothing.filter_taubin(
                mesh, lamb=0.05, nu=-0.53, iterations=10
            )
            verts = mesh.vertices
            faces = mesh.faces
        except Exception:
            pass

    cell = atoms.get_cell().array.T
    original_shape = np.asarray(dgrid.shape, dtype=float)
    shape_den = np.maximum(original_shape - 1.0, 1.0)
    verts_cart = np.dot(cell, ((verts - pad) / shape_den).T).T

    xx = np.arange(regions_padded.shape[0])
    yy = np.arange(regions_padded.shape[1])
    zz = np.arange(regions_padded.shape[2])
    interp = RGI(
        (xx, yy, zz),
        regions_padded,
        method="nearest",
        bounds_error=False,
        fill_value=0,
    )
    labels = interp(verts).astype(int)

    if select_labels is not None:
        select_set = set(int(x) for x in select_labels)
        vert_ok = np.isin(labels, list(select_set))
        face_mask = vert_ok[faces].all(axis=1)
        faces = faces[face_mask]

    if len(faces) == 0:
        return {
            "positions": [],
            "indices": [],
            "colors": [],
            "n_vertices": 0,
            "n_triangles": 0,
        }

    _DEFAULT_CSCALE = [
        [102, 197, 204], [246, 207, 113], [248, 156, 116],
        [220, 176, 242], [135, 197, 95], [158, 185, 243],
        [254, 136, 177], [201, 219, 116], [139, 224, 164],
        [180, 151, 231], [179, 179, 179],
    ]
    if cscale is None:
        cscale = _DEFAULT_CSCALE
    uniq = sorted([int(x) for x in np.unique(labels) if int(x) > 0])
    cmap = {rid: cscale[i % len(cscale)] for i, rid in enumerate(uniq)}
    colors = (
        np.asarray([cmap.get(int(r), [180, 180, 180]) for r in labels], dtype=float)
        / 255.0
    )

    return {
        "positions": verts_cart.reshape(-1).tolist(),
        "indices": faces.reshape(-1).tolist(),
        "colors": colors.reshape(-1).tolist(),
        "n_vertices": int(len(verts_cart)),
        "n_triangles": int(len(faces)),
    }


def _make_flat_dgrid_surface_trace(
    dgrid,
    atoms,
    isovalue=1.0,
    step_size=2,
    trace_id="surface_dgrid_flat",
    name="Dgrid surface (flat)",
    flat_color=None,
):
    """Build a flat-colored Three.js mesh trace from the raw distance grid.

    Parameters
    ----------
    dgrid : np.ndarray
        Distance grid.
    atoms : ase.Atoms
        Unit-cell reference.
    isovalue : float
        Marching-cubes isosurface level.
    step_size : int
        Marching-cubes step size.
    trace_id, name
        Trace identifiers.
    flat_color : list[float], optional
        RGB color (0-1).  Default is skyblue ``[0.53, 0.81, 0.92]``.

    Returns
    -------
    dict
        Three.js mesh payload.
    """
    verts, faces, _normals, _vals = marching_cubes(
        dgrid, level=isovalue, step_size=max(step_size, 1)
    )
    if len(faces) == 0:
        return {
            "positions": [],
            "indices": [],
            "colors": [],
            "n_vertices": 0,
            "n_triangles": 0,
        }

    cell = atoms.get_cell().array.T
    shape_den = np.maximum(np.asarray(dgrid.shape, dtype=float) - 1.0, 1.0)
    verts_cart = np.dot(cell, (verts / shape_den).T).T

    if flat_color is None:
        flat_color = [0.53, 0.81, 0.92]  # skyblue
    flat_color = np.asarray(flat_color, dtype=float)
    colors = np.tile(flat_color, (len(verts_cart), 1))

    return {
        "positions": verts_cart.reshape(-1).tolist(),
        "indices": faces.reshape(-1).tolist(),
        "colors": colors.reshape(-1).tolist(),
        "n_vertices": int(len(verts_cart)),
        "n_triangles": int(len(faces)),
    }


def _make_padding_maxima_payload(
    dgrid,
    padding_maxima,
    atoms,
    colors=None,
    cscale=None,
    opacity=0.7,
):
    """Build a Three.js points payload for padding maxima (outside central cell)."""
    if len(padding_maxima) == 0:
        return {"positions": [], "colors": [], "sizes": [], "n_points": 0}

    cell = atoms.get_cell().array.T
    shape_den = np.maximum(np.asarray(dgrid.shape, dtype=float) - 1.0, 1.0)
    frac = padding_maxima / shape_den
    pos = np.dot(cell, frac.T).T

    wrapped = np.mod(padding_maxima, np.asarray(dgrid.shape, dtype=int)).astype(int)
    radii = dgrid[tuple(wrapped.T)]
    sizes = np.clip(radii * 1.5, 1.5, 6.0).astype(float).tolist()

    if colors is not None:
        per_point_colors = np.asarray(colors, dtype=float)
    else:
        _DEFAULT_CSCALE = [
            [102, 197, 204], [246, 207, 113], [248, 156, 116],
            [220, 176, 242], [135, 197, 95], [158, 185, 243],
            [254, 136, 177], [201, 219, 116], [139, 224, 164],
            [180, 151, 231], [179, 179, 179],
        ]
        if cscale is None:
            cscale = _DEFAULT_CSCALE
        labels = dgrid[tuple(wrapped.T)].astype(int)
        uniq = sorted([int(x) for x in np.unique(labels) if int(x) > 0])
        cmap = {rid: cscale[i % len(cscale)] for i, rid in enumerate(uniq)}
        per_point_colors = (
            np.asarray([cmap.get(int(r), [180, 180, 180]) for r in labels], dtype=float)
            / 255.0
        )

    return {
        "positions": pos.reshape(-1).tolist(),
        "colors": per_point_colors.reshape(-1).tolist(),
        "sizes": sizes,
        "n_points": int(len(padding_maxima)),
    }


# ---------------------------------------------------------------------------
# Reusable builder
# ---------------------------------------------------------------------------

def build_padded_threejs_scene(
    cif_path,
    output_path=None,
    spacing=0.25,
    isovalue=1.0,
    camera_eye=(2.0, 0.0, 1.35),
    camera_center=(0.0, 0.0, 0.0),
    taubin_lambda=0.05,
    **seg_kwargs,
):
    """Build a Three.js HTML scene with padded transparent segmentation surfaces.

    Parameters
    ----------
    cif_path : str or Path
        Path to the CIF file.
    output_path : str or Path, optional
        Destination HTML path.  Defaults to ``<cif_stem>_padded_threejs.html``
        in the current working directory.
    spacing : float
        Distance-grid voxel spacing (Angstrom).  Default 0.25 for finer meshes.
    isovalue : float
        Marching-cubes isosurface level.
    camera_eye, camera_center
        Three.js camera position / look-at.
    taubin_lambda : float
        Taubin smoothing lambda parameter (default 0.05).
    **seg_kwargs
        Any keyword argument accepted by ``do_periodic_segmentation_rust``
        (e.g. ``mask_thickness``, ``h``, ``threshold_abs``, ``min_distance``,
        ``extra_padding``, ``preserve_threshold``, ``compactness``, …).
        The builder always forces ``return_padded_regions=True``,
        ``return_dgrid=True`` and ``compute_connections=True``.

    Returns
    -------
    Path
        Path to the generated HTML file.
    """
    cif_path = Path(cif_path)
    if not cif_path.exists():
        raise FileNotFoundError(f"CIF not found: {cif_path}")

    if output_path is None:
        output_path = Path.cwd() / f"{cif_path.stem}_padded_threejs.html"
    output_path = Path(output_path)

    print(f"Reading {cif_path} ...")
    atoms = read(str(cif_path))
    atoms_updated = mgr.update_with_covalent_bonds(atoms)

    print(f"Computing distance grid (Rust, spacing={spacing}) ...")
    dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=spacing)

    print(f"Running periodic segmentation (Rust) with seg_kwargs={seg_kwargs} ...")
    results = mgr.do_periodic_segmentation_rust(
        dgrid,
        atoms,
        spacing=spacing,
        return_padded_regions=True,
        return_dgrid=True,
        compute_connections=True,
        **seg_kwargs,
    )

    regions = results["regions"]
    regions_pbc = results["regions_pbc"]
    regions_padded = results["regions_padded"]
    all_maxima = results["all_maxima"]
    maxima_pbc = results["maxima_pbc"]
    connections = results["connections"]
    unique_regions = results["unique_regions"]

    _extra_padding = seg_kwargs.get("extra_padding", 2.0)
    padding_thickness = int(
        mgr.calculate_padding_thickness(dgrid, atoms, spacing)
        + np.ceil(_extra_padding / spacing)
    )
    print(f"Padding thickness: {padding_thickness} voxels")

    print("Building MovieFigure ...")
    fig = MovieFigure()
    fig.set_camera(eye=camera_eye, center=camera_center)

    # 1) Framework
    fig.add_framework(atoms_updated)
    # 1b) Aromatic rings
    fig.add_aromatic_rings(atoms_updated)
    # 2) Cell
    fig.add_cell(atoms)

    # 3) Standard regions surface
    fig.add_surface(dgrid, regions, atoms, isovalue=isovalue, step_size=2)

    # 4) Regions PBC surface
    trace_regions_pbc = make_surface_trace(
        dgrid,
        regions_pbc,
        atoms,
        isovalue=isovalue,
        step_size=2,
        trace_id="surface_regions_pbc",
        name="Regions (PBC)",
    )
    fig.add_trace(trace_regions_pbc)

    # 5) Flat skyblue dgrid surface (no segmentation colors — single flat color)
    payload_dgrid_flat = _make_flat_dgrid_surface_trace(
        dgrid,
        atoms,
        isovalue=isovalue,
        step_size=2,
        trace_id="surface_dgrid_flat",
        name="Dgrid surface (flat)",
    )
    trace_dgrid_flat = _make_trace(
        trace_id="surface_dgrid_flat",
        trace_type="mesh",
        name="Dgrid surface (flat)",
        data=payload_dgrid_flat,
        style={"opacity": 0.6, "side": "DoubleSide"},
    )
    fig.add_trace(trace_dgrid_flat)
    fig.groups["dgrid_flat_surface"] = {
        "traces": ["surface_dgrid_flat"],
        "label": "Dgrid surface (flat)",
        "default_off": True,
    }

    # 6) Pore graph (maxima + connections + windows)
    seg = {
        "dgrid": dgrid,
        "regions": regions,
        "all_maxima": all_maxima,
        "maxima_pbc": maxima_pbc,
        "connections": connections,
    }
    fig.add_maxima(seg, atoms, hide_periodic=True)
    fig.add_connections(seg, atoms, hide_periodic=True)
    fig.add_windows(seg, hide_periodic=True)

    # 6b) Legend groups for pore-graph traces (so they appear as toggleable entries)
    fig.groups["maxima_internal"] = {
        "traces": ["maxima_internal"],
        "label": "Maxima (internal)",
        "default_off": False,
    }
    fig.groups["maxima_outside"] = {
        "traces": ["maxima_outside"],
        "label": "Maxima (periodic)",
        "default_off": True,
    }
    fig.groups["maxima_pbc"] = {
        "traces": ["maxima_pbc"],
        "label": "Maxima (PBC reps)",
        "default_off": False,
    }
    fig.groups["connections_internal"] = {
        "traces": ["connections_internal", "connections_internal_nodes"],
        "label": "Connections (internal)",
        "default_off": False,
    }
    fig.groups["connections_periodic"] = {
        "traces": ["connections_periodic", "connections_periodic_nodes"],
        "label": "Connections (periodic)",
        "default_off": True,
    }
    fig.groups["windows_internal"] = {
        "traces": ["windows_internal"],
        "label": "Windows (internal)",
        "default_off": False,
    }
    fig.groups["windows_periodic"] = {
        "traces": ["windows_periodic"],
        "label": "Windows (periodic)",
        "default_off": True,
    }

    # 7) PADDED REGIONS (non-PBC) — only the padding envelope, transparent
    all_padded_labels = set(int(x) for x in np.unique(regions_padded) if x > 0)
    central_labels = set(int(x) for x in np.unique(regions) if x > 0)
    labels_in_padding = sorted(all_padded_labels - central_labels)
    print(f"Padding-only labels: {len(labels_in_padding)}")

    payload_padded_only = _make_padded_surface_trace(
        dgrid,
        regions_padded,
        atoms,
        padding_thickness=padding_thickness,
        isovalue=isovalue,
        step_size=2,
        opacity=0.2,
        trace_id="surface_padded_only",
        name="Padded regions (only padding)",
        select_labels=labels_in_padding,
    )
    trace_padded_only = _make_trace(
        trace_id="surface_padded_only",
        trace_type="mesh",
        name="Padded regions (only padding)",
        data=payload_padded_only,
        style={"opacity": 0.2, "side": "DoubleSide"},
    )
    fig.add_trace(trace_padded_only)

    # 8) PADDED REGIONS (PBC-aware) — wrap regions_pbc into padding, show all, transparent
    padded_regions_pbc = np.pad(regions_pbc, padding_thickness, mode="wrap")
    payload_padded_pbc = _make_padded_surface_trace(
        dgrid,
        padded_regions_pbc,
        atoms,
        padding_thickness=padding_thickness,
        isovalue=isovalue,
        step_size=2,
        opacity=0.2,
        trace_id="surface_padded_pbc",
        name="Padded regions (PBC)",
        select_labels=None,
    )
    trace_padded_pbc = _make_trace(
        trace_id="surface_padded_pbc",
        trace_type="mesh",
        name="Padded regions (PBC)",
        data=payload_padded_pbc,
        style={"opacity": 0.2, "side": "DoubleSide"},
    )
    fig.add_trace(trace_padded_pbc)

    # 9) GROUPS / LEGEND entries so padded surfaces are toggleable
    fig.groups["padded_regions_only"] = {
        "traces": ["surface_padded_only"],
        "label": "Padded regions (only padding)",
        "default_off": True,
    }
    fig.groups["padded_regions_pbc"] = {
        "traces": ["surface_padded_pbc"],
        "label": "Padded regions (PBC)",
        "default_off": True,
    }

    # 10) PADDING MAXIMA
    maxima_info = results["maxima_info"]
    maxima_pbc_dict = results.get("maxima_pbc_dict", {})
    n_central_regions = int(unique_regions.max())
    padding_maxima = all_maxima[n_central_regions - 1 :]
    mi_padding = maxima_info[n_central_regions - 1 :] if maxima_info is not None else None

    _ALPHABET_CSCALE = [
        [102, 197, 204], [246, 207, 113], [248, 156, 116],
        [220, 176, 242], [135, 197, 95], [158, 185, 243],
        [254, 136, 177], [201, 219, 116], [139, 224, 164],
        [180, 151, 231], [179, 179, 179],
    ]

    if len(padding_maxima) > 0 and mi_padding is not None:
        # 9a) Distinct region-label colours for padding maxima
        region_labels = mi_padding[:, 0].astype(int)
        uniq = sorted([int(x) for x in np.unique(region_labels) if int(x) > 0])
        cmap = {rid: _ALPHABET_CSCALE[i % len(_ALPHABET_CSCALE)] for i, rid in enumerate(uniq)}
        pm_colors = (
            np.asarray([cmap.get(int(r), [180, 180, 180]) for r in region_labels], dtype=float)
            / 255.0
        )
        _pm_payload = _make_padding_maxima_payload(
            dgrid, padding_maxima, atoms, colors=pm_colors, opacity=0.7
        )
        trace_padding_maxima = _make_trace(
            trace_id="padding_maxima",
            trace_type="points",
            name="Padding maxima",
            data=_pm_payload,
            style={"size": 2.5},
            visible=True,
        )
        fig.add_trace(trace_padding_maxima)
        fig.groups["padding_maxima"] = {
            "traces": ["padding_maxima"],
            "label": "Padding maxima",
            "default_off": True,
        }

        # 9b) PBC-group colours for padding maxima
        pbc_labels = []
        for m in mi_padding[:, 1:4]:
            key = (int(m[0]), int(m[1]), int(m[2]))
            lbl = maxima_pbc_dict.get(key)
            if lbl is None:
                lbl = maxima_pbc_dict.get(f"{key[0]},{key[1]},{key[2]}")
            if lbl is None:
                lbl = maxima_pbc_dict.get(str(key))
            pbc_labels.append(int(lbl) if lbl is not None else 0)

        pbc_labels = np.asarray(pbc_labels, dtype=int)
        if np.any(pbc_labels > 0):
            uniq_pbc = sorted([int(x) for x in np.unique(pbc_labels) if int(x) > 0])
            cmap_pbc = {rid: _ALPHABET_CSCALE[i % len(_ALPHABET_CSCALE)] for i, rid in enumerate(uniq_pbc)}
            pm_pbc_colors = (
                np.asarray([cmap_pbc.get(int(r), [180, 180, 180]) for r in pbc_labels], dtype=float)
                / 255.0
            )
            _pm_pbc_payload = _make_padding_maxima_payload(
                dgrid, padding_maxima, atoms, colors=pm_pbc_colors, opacity=0.7
            )
            trace_padding_maxima_pbc = _make_trace(
                trace_id="padding_maxima_pbc",
                trace_type="points",
                name="Padding maxima (PBC colours)",
                data=_pm_pbc_payload,
                style={"size": 2.5},
                visible=True,
            )
            fig.add_trace(trace_padding_maxima_pbc)
            fig.groups["padding_maxima_pbc"] = {
                "traces": ["padding_maxima_pbc"],
                "label": "Padding maxima (PBC colours)",
                "default_off": True,
            }

    # --- Export ---
    print(f"Writing HTML to {output_path} ...")
    fig.to_html(str(output_path))
    print(f"Done.  ({output_path.stat().st_size / 1e6:.1f} MB)")
    return output_path


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Build a Three.js HTML scene with padded transparent segmentation surfaces."
    )
    parser.add_argument(
        "--cif",
        required=True,
        help="Path to the input CIF file (e.g. examples/TestCIFs/P1-IRMOF-1.cif)",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Output HTML path (default: <cif_stem>_padded_threejs.html)",
    )
    parser.add_argument(
        "--spacing",
        type=float,
        default=0.25,
        help="Distance-grid spacing in Angstrom (default: 0.25)",
    )
    parser.add_argument(
        "--isovalue",
        type=float,
        default=1.0,
        help="Marching-cubes isosurface level (default: 1.0)",
    )
    parser.add_argument(
        "--taubin-lambda",
        type=float,
        default=0.05,
        help="Taubin smoothing lambda (default: 0.05)",
    )
    # --- do_periodic_segmentation_rust kwargs ---
    parser.add_argument(
        "--mask-thickness",
        type=float,
        default=0.5,
        help="Void-space mask threshold (default: 0.5)",
    )
    parser.add_argument(
        "--h",
        type=float,
        default=0.2,
        help="H-minima transform height for merging maxima (default: 0.2)",
    )
    parser.add_argument(
        "--threshold-abs",
        type=float,
        default=1.5,
        help="Absolute threshold for maxima detection (default: 1.5)",
    )
    parser.add_argument(
        "--min-distance",
        type=float,
        default=1.5,
        help="Minimum separation between maxima / windows (default: 1.5)",
    )
    parser.add_argument(
        "--extra-padding",
        type=float,
        default=2.0,
        help="Extra padding thickness in Angstrom (default: 2.0)",
    )
    parser.add_argument(
        "--preserve-threshold",
        type=float,
        default=3.0,
        help="Preserve-threshold for maxima merging (default: 3.0)",
    )
    parser.add_argument(
        "--smooth-threshold",
        type=float,
        default=6.0,
        help="Smooth-threshold for maxima merging (default: 6.0)",
    )
    parser.add_argument(
        "--max-smooth-radius",
        type=float,
        default=2.0,
        help="Max smooth radius for maxima merging (default: 2.0)",
    )
    parser.add_argument(
        "--wall-thickness",
        type=float,
        default=0.05,
        help="Wall-thickness for wall-maxima removal (default: 0.05)",
    )
    parser.add_argument(
        "--compactness",
        type=float,
        default=0.0,
        help="Watershed compactness (default: 0.0)",
    )
    parser.add_argument(
        "--min-region-size",
        type=float,
        default=0.001,
        help="Minimum region size fraction (default: 0.001)",
    )
    parser.add_argument(
        "--watershed-backend",
        type=str,
        default="rust_bucket",
        help="Watershed backend: rust_bucket | rust_det | rust_det_parallel | skimage (default: rust_bucket)",
    )
    parser.add_argument(
        "--watershed-connectivity",
        type=int,
        default=2,
        help="Watershed connectivity (1=6-NN, 2=26-NN) (default: 2)",
    )
    parser.add_argument(
        "--minimum-window-separation",
        type=float,
        default=None,
        help="Minimum window separation (default: spacing*3)",
    )
    parser.add_argument(
        "--disable-tqdm",
        action="store_true",
        help="Disable tqdm progress bars",
    )

    args = parser.parse_args()

    seg_kwargs = {
        "mask_thickness": args.mask_thickness,
        "h": args.h,
        "threshold_abs": args.threshold_abs,
        "min_distance": args.min_distance,
        "extra_padding": args.extra_padding,
        "preserve_threshold": args.preserve_threshold,
        "smooth_threshold": args.smooth_threshold,
        "max_smooth_radius": args.max_smooth_radius,
        "wall_thickness": args.wall_thickness,
        "compactness": args.compactness,
        "min_region_size": args.min_region_size,
        "watershed_backend": args.watershed_backend,
        "watershed_connectivity": args.watershed_connectivity,
        "disable_tqdm": args.disable_tqdm,
    }
    if args.minimum_window_separation is not None:
        seg_kwargs["minimum_window_separation"] = args.minimum_window_separation

    build_padded_threejs_scene(
        cif_path=args.cif,
        output_path=args.output,
        spacing=args.spacing,
        isovalue=args.isovalue,
        taubin_lambda=args.taubin_lambda,
        **seg_kwargs,
    )


if __name__ == "__main__":
    main()
