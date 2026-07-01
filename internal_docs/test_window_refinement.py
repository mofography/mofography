#%%
"""Ad-hoc script to compare window refinement prototypes."""

import logging
from pathlib import Path
from typing import Iterable, List

import numpy as np
import plotly.graph_objects as go
from plotly.offline import plot as plotly_plot
from ase.io import read
import plotly.express as px
from mofography.window_refinement import (
    WindowRefinementOptions,
    refine_windows_via_morph_opening,
)

import mofography as mgr

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

TARGET_CIF = "NU-1000_S1.cif"


def build_segmentation_result(path_to_cif: Path):
    mat_atoms = read(path_to_cif)
    spacing = 0.5
    params = dict(
        mask_thickness=0.5,
        h=0.2,
        threshold_abs=1.5,
        min_distance=1.5,
        minimum_window_separation=1.5,
        wall_thickness=0.08,
        extra_padding=10,
        return_dgrid=True,
        return_padded_regions=True,
        compute_connections=True,
    )
    logger.info("Building distance grid (Rust backend)")
    dgrid = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=spacing)
    logger.info("Running periodic segmentation (includes connections)")
    results = mgr.do_periodic_segmentation_rust(
        dgrid,
        mat_atoms,
        spacing=spacing,
        **params,
    )
    results["cell_matrix"] = mat_atoms.cell.array.T
    return results, mat_atoms


def _minimal_image_distance(a, b, cell_matrix, cell_inv):
    frac_a = cell_inv @ a
    frac_b = cell_inv @ b
    dfrac = frac_b - frac_a
    dfrac -= np.round(dfrac)
    delta = cell_matrix @ dfrac
    return float(np.linalg.norm(delta))


def compute_deltas(baseline, refined, cell_matrix):
    center_deltas = []
    radius_deltas = []
    cell_inv = np.linalg.inv(cell_matrix)
    for base_conn, new_conn in zip(baseline, refined):
        base_centers = np.asarray(base_conn[3])
        new_centers = np.asarray(new_conn[3])
        base_radii = np.asarray(base_conn[4])
        new_radii = np.asarray(new_conn[4])
        assert (
            base_centers.shape == new_centers.shape
        ), "Connection window counts mismatch"
        for bc, nc in zip(base_centers, new_centers):
            center_deltas.append(
                _minimal_image_distance(bc, nc, cell_matrix, cell_inv)
            )
        for br, nr in zip(base_radii, new_radii):
            radius_deltas.append(float(abs(br - nr)))
    return np.array(center_deltas), np.array(radius_deltas)


def flatten_radii(connections: Iterable[np.ndarray]) -> np.ndarray:
    values: List[float] = []
    for conn in connections:
        values.extend(np.asarray(conn[4], dtype=float).tolist())
    return np.asarray(values, dtype=float)


def plot_histogram(
    deltas: np.ndarray,
    *,
    title: str,
    xlabel: str,
    outfile: Path,
):
    fig = go.Figure()
    if deltas.size:
        fig.add_trace(
            go.Histogram(
                x=deltas,
                nbinsx=40,
                name="Morph opening",
                opacity=0.75,
            )
        )
    fig.update_layout(
        title=title,
        xaxis_title=xlabel,
        yaxis_title="Count",
    )
    plotly_plot(
        fig,
        filename=str(outfile),
        auto_open=False,
        include_plotlyjs="cdn",
    )
    logger.info(f"Saved plot to {outfile}")


def plot_parity(
    baseline_r: np.ndarray,
    refined_r: np.ndarray,
    *,
    title: str,
    outfile: Path,
):
    fig = go.Figure(
        go.Scatter(
            x=baseline_r,
            y=refined_r,
            mode="markers",
            marker=dict(size=6, opacity=0.7),
            name="windows",
        )
    )
    if baseline_r.size and refined_r.size:
        line_min = min(baseline_r.min(), refined_r.min())
        line_max = max(baseline_r.max(), refined_r.max())
    else:
        line_min, line_max = 0.0, 1.0
    fig.add_trace(
        go.Scatter(
            x=[line_min, line_max],
            y=[line_min, line_max],
            mode="lines",
            line=dict(color="black", dash="dash"),
            name="y = x",
        )
    )
    fig.update_layout(
        title=title,
        xaxis_title="Baseline radius (Å)",
        yaxis_title="Refined radius (Å)",
    )
    plotly_plot(
        fig,
        filename=str(outfile),
        auto_open=False,
        include_plotlyjs="cdn",
    )
    logger.info(f"Saved parity plot to {outfile}")


def summarize_deltas(label: str, center_delta: np.ndarray, radius_delta: np.ndarray):
    summary = {
        f"{label}_center_mean": float(center_delta.mean()) if center_delta.size else 0.0,
        f"{label}_center_min": float(center_delta.min(initial=0.0)),
        f"{label}_center_max": float(center_delta.max(initial=0.0)),
        f"{label}_center_median": float(np.median(center_delta))
        if center_delta.size
        else 0.0,
        f"{label}_radius_mean": float(radius_delta.mean()) if radius_delta.size else 0.0,
        f"{label}_radius_min": float(radius_delta.min(initial=0.0)),
        f"{label}_radius_max": float(radius_delta.max(initial=0.0)),
        f"{label}_radius_median": float(np.median(radius_delta))
        if radius_delta.size
        else 0.0,
    }
    return summary


def trim_connections(connections, unique_regions):
    return [
        c for c in connections if np.all(np.isin(c[0], unique_regions))
    ]


def scene_plot(
    results: dict,
    mat_atoms,
    out_path: Path,
    *,
    connections_key: str = "connections",
    isovals=(1.0, 1.5),
):
    dgrid = results["dgrid"]
    regions = results["regions"]
    regions_pbc = results["regions_pbc"]
    maxima_info = results["maxima_info"]
    pbc_groups = results["pbc_groups"]
    connections = results[connections_key]
    all_maxima = results["all_maxima"]
    all_maxima_radii = dgrid[tuple(maxima_info[:, 1:4].T)]
    maxima_pbc = results["maxima_pbc"]
    unique_regions = results["unique_regions"]
    maxima = all_maxima[: regions.max()]
    maxima_radii = all_maxima_radii[: regions.max()]
    connections_trimmed = trim_connections(connections, unique_regions)

    cell_trace = mgr.plot_cell(mat_atoms, color="goldenrod")
    cscale = px.colors.qualitative.Pastel
    isosurface_segmentation = mgr.plot_surface_dgrid(
        dgrid,
        stylize=True,
        isovalue=isovals[0],
        region_labels=regions,
        ase_atoms=mat_atoms,
        name="Regions",
        opacity=0.9,
        cscale=cscale,
    )
    isosurface_pbc = mgr.plot_surface_dgrid(
        dgrid,
        stylize=True,
        isovalue=isovals[0],
        region_labels=regions_pbc,
        ase_atoms=mat_atoms,
        name="Regions PBC",
        opacity=0.4,
        cscale=cscale,
    )
    maxima_trace = mgr.plot_maxima(
        dgrid,
        maxima,
        mat_atoms,
        marker_size=10,
        show_spheres=False,
        marker_colorscale=cscale,
        name="Maxima (distinct)",
    )
    maxima_trace_pbc = mgr.plot_maxima(
        dgrid,
        maxima_pbc,
        mat_atoms,
        show_spheres=False,
        marker_size=7.5,
        marker_colorscale=cscale,
        name="Maxima (pbc)",
    )
    connection_traces = mgr.plot_connections_as_it_is(
        connections_trimmed,
        all_maxima,
        regions.shape,
        mat_atoms,
        cscale=cscale,
        plot_disconnected=False,
        plot_periodic=True,
    )
    windows_trace = mgr.plot_windows(
        connections_trimmed,
        name="Windows",
        marker_color="yellow",
        marker_size=8,
        marker_opacity=0.8,
    )

    fig = go.Figure(
        data=[
            cell_trace,
            *connection_traces,
            isosurface_segmentation,
            isosurface_pbc,
            maxima_trace,
            maxima_trace_pbc,
            windows_trace,
        ],
        layout=mgr.white_layout,
    )
    plotly_plot(
        fig,
        filename=str(out_path),
        auto_open=False,
        include_plotlyjs="cdn",
    )
    logger.info("Saved scene to %s", out_path)


if __name__ == "__main__":
    path_to_cif = (
        Path("~/repos/mofography/examples/TestCIFs").expanduser() / TARGET_CIF
    )
    assert path_to_cif.exists(), f"CIF not found: {path_to_cif}"
    results, mat_atoms = build_segmentation_result(path_to_cif)
    baseline_connections = results["connections"]

    options = WindowRefinementOptions(search_radius_vox=8)
    refined_morph = refine_windows_via_morph_opening(
        results,
        cell_matrix=results["cell_matrix"],
        options=options,
    )

    cell_matrix = results["cell_matrix"]
    morph_shift, morph_radius = compute_deltas(
        baseline_connections, refined_morph, cell_matrix
    )

    summary = summarize_deltas("morph", morph_shift, morph_radius)
    print("Summary:", summary)

    plot_dir = Path(__file__).with_name("window_refinement_plots")
    plot_dir.mkdir(exist_ok=True)

    plot_histogram(
        morph_shift,
        title="Window center shifts vs. baseline (Å)",
        xlabel="‖Δcenter‖ (Å)",
        outfile=plot_dir / "center_shift_hist.html",
    )
    plot_histogram(
        morph_radius,
        title="Window radius changes vs. baseline (Å)",
        xlabel="|Δradius| (Å)",
        outfile=plot_dir / "radius_delta_hist.html",
    )

    base_radii = flatten_radii(baseline_connections)
    morph_radii = flatten_radii(refined_morph)

    plot_parity(
        base_radii,
        morph_radii,
        title="Baseline vs. morph-opening radii",
        outfile=plot_dir / "parity_morph.html",
    )

    results["connections_morph_opening"] = refined_morph

    scene_plot(
        results,
        mat_atoms,
        plot_dir / "scene_baseline.html",
        connections_key="connections",
    )
    scene_plot(
        results,
        mat_atoms,
        plot_dir / "scene_morph_opening.html",
        connections_key="connections_morph_opening",
    )