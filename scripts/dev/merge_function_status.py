"""Merge per-module status proposals + owner overrides into docs/function_status.yml.

Proposals: internal_docs/docs_redesign/status/*.yml (name -> status/note/use/backend).
Output keys are "<package>.<name>" so same-named functions in Plotly and Three.js modules stay distinct.
Run from repo root:  python scripts/dev/merge_function_status.py
"""
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
SDIR = ROOT / "internal_docs" / "docs_redesign" / "status"
OUT = ROOT / "docs" / "function_status.yml"

FILES = {  # package -> proposal files
    "mofography.compute": ["compute_a_f", "compute_g_z"],
    "mofography.visualize": ["visualize"], "mofography.viz_backend": ["visualize"],
    "mofography.analyze": ["analyze_data_helper"], "mofography.data": ["analyze_data_helper"],
    "mofography.helper": ["analyze_data_helper"], "mofography.vdw_radii": ["analyze_data_helper"],
    "mofography.cif_io": ["analyze_data_helper"],
    "mofography.rust_fun": ["rust_workflow_typing"], "mofography.workflow": ["rust_workflow_typing"],
    "mofography.pocket_typing": ["rust_workflow_typing"],
    "mofography.pocket_framework_env": ["rust_workflow_typing"],
    "mofography.widest_path_resegment": ["rust_workflow_typing"],
    "mofography.window_refinement": ["rust_workflow_typing"],
    "mofography.threejs_viz": ["threejs_movies"], "mofography.threejs_plotting": ["threejs_movies"],
    "mofography.movies": ["threejs_movies"], "mofography.export_mp4": ["threejs_movies"],
    "mofography.winding_percolation": ["winding_raspa"], "mofography.percolation_results": ["winding_raspa"],
    "mofography.raspa_helper": ["winding_raspa"],
}

R, A, E, L, D, U = "recommended", "alternate", "experimental", "legacy", "deprecated", "utility"
W = "find_pore_limiting_radius_winding"
# Owner decisions (name -> (status, note, use)); applied to every module that has the name.
OVERRIDES = {
    "add_clustering_to_rag": (L, "Older CSD-based clustering; env_cluster typing is more accurate.", "add_env_cluster_pocket_types_to_rag"),
    "add_clustering_to_rag_energy": (L, "Older energy-based clustering; env_cluster typing is the default.", "add_env_cluster_pocket_types_to_rag"),
    "add_clustering_to_rag_energy_js": (L, "Older energy+JS clustering, less accurate than env_cluster.", "add_env_cluster_pocket_types_to_rag"),
    "add_clustering_to_rag_geometry_js": (L, "Older geometry+JS clustering, less accurate than env_cluster.", "add_env_cluster_pocket_types_to_rag"),
    "add_clustering_to_rag_pca": (D, "PCA-based clustering is no longer used.", "add_env_cluster_pocket_types_to_rag"),
    "add_csd_to_rag_with_pbc": (L, "CSD descriptor is old and simple; env_cluster typing supersedes it.", "add_env_cluster_pocket_types_to_rag"),
    "add_pbc_to_rag": (U, "Harmless helper that adds periodic-image info to a RAG; rarely needed.", None),
    "add_env_cluster_pocket_types_to_rag": (R, "Recommended default for pocket typing on a RAG.", None),
    "compute_percolation_history": (A, "Lattice method, kept for comparison; winding is the default.", "compute_percolation_history_winding"),
    "compute_percolation_history2": (A, "Lattice method, kept for comparison; winding is the default.", "compute_percolation_history_winding"),
    "compute_percolation_history_from_results": (A, "Lattice method, kept for comparison; winding is the default.", "compute_percolation_history_winding"),
    "compute_percolation_history_from_results2": (A, "Lattice method, kept for comparison; winding is the default.", "compute_percolation_history_winding"),
    "compute_percolation_history_from_results_lattice": (A, "Lattice method, kept for comparison; winding is the default.", "compute_percolation_history_winding"),
    "do_periodic_segmentation": (D, "Slow Python segmentation; not kept up to date.", "do_periodic_segmentation_rust"),
    "do_periodic_segmentation_basic": (D, "Slow Python variant, not kept up to date.", "do_periodic_segmentation_rust"),
    "find_isomorphic_path_groups": (E, "Experimental path grouping.", None),
    "find_paths": (E, "Experimental path search.", None),
    "find_unique_sgraphs": (E, "Experimental structure-graph deduplication.", None),
    "find_pore_limiting_radius": (L, "Lattice method, not recommended.", W),
    "find_pore_limiting_radius2": (L, "Lattice method, not recommended.", W),
    "get_energy_interpolator_tricubic": (R, "Most accurate energy interpolator; used for configurational-integral screening in the paper.", None),
    "get_pld_for_rag": (A, "RAG-based PLD; the winding functions are the default.", W),
    "get_pld_from_percolated_paths_path_windows": (D, "Old RAG-path PLD helper; use the winding PLD instead.", W),
    "get_pore_dia_mgr": (A, "Lattice variant, not the default.", "get_pore_dia_mgr_winding"),
    "get_pore_dia_mgr2": (A, "Lattice variant, not the default.", "get_pore_dia_mgr_winding"),
    "get_wall_windows2": (U, "Wall-window detection used inside get_connections_for_rag.", None),
    "get_wall_windows3": (D, "Unfinished stub that does not run.", "get_wall_windows2"),
    "plot_pore_graph": (L, "Old pore-graph plot for node_coords graphs.", "plot_connections_as_it_is"),
    "plot_rag": (A, "Optional alternate RAG plot.", "plot_connections_as_it_is"),
    "plot_connections_as_it_is": (R, "Go-to connection/RAG plot.", None),
    "plot_connections_as_it_is_3js": (R, "Go-to Three.js connection plot.", None),
    "rag_from_connections": (R, "Go-to RAG builder.", None),
}
OVERRIDE_PKGS_EXCLUDE = {  # name -> packages where override must NOT apply
    "plot_connections_as_it_is": {"mofography.threejs_viz"},
}


def main():
    cfg = yaml.safe_load((ROOT / "docs" / "_quarto.yml").read_text())["quartodoc"]
    props = {p.stem: yaml.safe_load(p.read_text()) or {} for p in SDIR.glob("*.yml")}
    out, missing = {}, []
    for s in cfg["sections"]:
        pkg = s["package"]
        for n in s["contents"]:
            n = n if isinstance(n, str) else n["name"]
            e = next((props[f][n] for f in FILES.get(pkg, []) if n in props.get(f, {})), None)
            if e is None:
                missing.append(f"{pkg}.{n}")
                e = {"status": U, "note": ""}
            rec = {"status": e["status"], "note": (e.get("note") or "").strip()}
            if e.get("use"):
                rec["use"] = e["use"]
            if n in OVERRIDES and pkg not in OVERRIDE_PKGS_EXCLUDE.get(n, ()):
                st, note, use = OVERRIDES[n]
                rec = {"status": st, "note": note}
                if use:
                    rec["use"] = use
            out[f"{pkg}.{n}"] = rec
    OUT.write_text("# Generated by scripts/dev/merge_function_status.py; edit OVERRIDES there.\n"
                   + yaml.safe_dump(out, sort_keys=True, allow_unicode=True, width=110))
    from collections import Counter
    print(Counter(v["status"] for v in out.values()), "\nmissing:", missing)


if __name__ == "__main__":
    main()
