# Flexible framework movies (Three.js HTML)

This note describes the **implemented** path in `mofography/movies.py` for exporting an interactive HTML viewer when the **framework coordinates change each frame** (flexible RASPA-style trajectory) and you want **matched adsorbate / component snapshots** on the same frame indices.

Design background lives in `docs/movie_viewer_module_plan.md`; this file is a **how-to** for the current API.

## When to use which entry point

| Goal | Entry point |
|------|-------------|
| Fixed framework from a CIF + trajectory of guests only | `make_movie_from_trajectory(...)` or `build_fixed_framework_movie(...)` |
| Flexible framework + segmentation/pore graph **per sparse frame** + one or more guest/component movies aligned by index | `make_sparse_multicomponent_movie(...)` |

## Requirements

1. **`framework_movie_path`** — multi-frame RASPA PDB (or ASE-readable path with `movie_format="ase"`) whose models are **only the framework** atoms, **same frame count** as the component movies you pass.
2. **`component_movie_paths`** — dict mapping a short label → path. Each movie must have the **same number of frames** as the framework movie (after any shared `frame_stride` / index selection). Typical RASPA outputs: `Movie_*_allcomponents.pdb` or per-component `Movie_*_component_X_*.pdb` files from the **same** `Run_*` folder.
3. **`cif_path`** — used for metadata / consistency checks in the builder; the **geometry** in the HTML comes from the **framework movie** frames when `framework_movie_path` is set.
4. **`snapshot_indices`** — optional explicit list of **source** frame indices to keep (sparse export). Segmentation (`dgrid`, regions, connections, surface mesh, maxima, windows) is computed **for each** selected index.

## Minimal example (single component, flexible framework)

```python
from pathlib import Path
import mofography as mgr

base = Path("outputs/github_raspa_test")  # or your data directory

out_html, payload = mgr.make_sparse_multicomponent_movie(
    component_movie_paths={
        "allcomponents": base / "Movie_IRMOF1_allcomponents.pdb",
    },
    framework_movie_path=base / "Movie_IRMOF1_frameworks.pdb",
    cif_path=base / "IRMOF-1.cif",
    out_html=base / "flexible_single_component_ir_mof1.html",
    snapshot_indices=[0, 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55],
    max_frames=15,
    movie_format="raspa_pdb",
    spacing=0.25,
    isovalue=1.0,
    surf_step_size=2,
)
print(out_html)
# payload["segmentation_series"] — one dict per selected frame (dgrid, regions, …)
```

Returned `payload` includes:

- `snapshot_indices` — indices actually embedded
- `segmentation_series` — list of `{ "frame_index", "segmentation", "summary" }` for analysis or saving `.npz` sidecars

## Public RASPA test assets (IRMOF-1 minimization)

From the **numat/RASPA2** mirror (same example folder):

- Framework movie:  
  `examples/Non-Basic/MinimizationIRMOF-1-NVT/Movies/System_0/Run_0/Movie_IRMOF-1_1.1.1_298.000000_0.000000_frameworks.pdb`
- All components:  
  `examples/Non-Basic/MinimizationIRMOF-1-NVT/Movies/System_0/Run_0/Movie_IRMOF-1_1.1.1_298.000000_0.000000_allcomponents.pdb`
- CIF: `structures/mofs/cif/IRMOF-1.cif`

Raw URLs (for `curl` or CI):

- https://raw.githubusercontent.com/numat/RASPA2/master/examples/Non-Basic/MinimizationIRMOF-1-NVT/Movies/System_0/Run_0/Movie_IRMOF-1_1.1.1_298.000000_0.000000_frameworks.pdb  
- https://raw.githubusercontent.com/numat/RASPA2/master/examples/Non-Basic/MinimizationIRMOF-1-NVT/Movies/System_0/Run_0/Movie_IRMOF-1_1.1.1_298.000000_0.000000_allcomponents.pdb  
- https://raw.githubusercontent.com/numat/RASPA2/master/structures/mofs/cif/IRMOF-1.cif  

## Trace styling (parity with `md_pdb_threejs_viewer.py`)

The Three.js template in `mofography/movies.py` reads **`trace.style` as a resolved object** on each trace (not only string keys into `DATA.styles`). Defaults match the devtools viewer:

- **Internal connections:** red (`0xFF0000`), solid, linewidth 16  
- **Periodic connections:** red, **dashed** (`dashSize` / `gapSize` like the devtools script), linewidth 16  
- **Windows (internal + periodic):** yellow `0xFFE600`, larger point size  
- **Internal maxima:** colored by **segmentation pocket / region** (same palette as the surface), marker size from cage radii (`6 + 5 * norm`) and capped by the viewer’s point size logic  

Guest / component spheres use distinct hues via `guests_<label>` trace ids.

## Related scripts

- **`devtools/scripts/md_pdb_threejs_viewer.py`** — legacy single-file payload (`framework`, `surface`, `movie`, …) and its own inline HTML. Prefer **`MovieFigure` / `make_sparse_multicomponent_movie`** for flexible + multi-trace workflows so logic stays in `movies.py`.

## See also

- [Aromatic ring detection](aromatic_ring_detection.qmd) — bond-length tolerances for ring highlights (Plotly & Three.js)
- `docs/movie_viewer_module_plan.md` — architecture and roadmap  
- `mofography/movies.py` — `read_movie_snapshots`, `compute_sparse_segmentation_series`, `make_sparse_multicomponent_movie`, `DEFAULT_STYLES`
