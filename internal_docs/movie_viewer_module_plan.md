# Movie Viewer Module Plan (Revised)

## 1) Goals and Scope

- Build a reusable Mofography movie module with Plotly-like ergonomics (`MovieFigure` + dict-based traces + layout/camera updates).
- Reuse existing `mofography/compute.py` and `mofography/analyze.py` logic for geometry/segmentation computation.
- Support two primary workflows:
  - **Fixed framework:** framework + pore traces remain static, adsorbates animate by frame.
  - **Flexible framework:** framework/surfaces/graph recomputed for all or selected snapshot indices.
- Separate heavy computation from visualization/export so preprocessing can run independently (and later in parallel).
- **Default bond handling:** use the existing `mofography.analyze.update_with_covalent_bonds` (PyMagin Voronoi + covalent radii tolerance) as the `auto` mode. Users may optionally pass precomputed bonds via `atoms.info['covalent_bonds']` or disable bond generation (`off`).

## 2) Proposed Architecture

Create `mofography/movie/` with clear layers:

- **I/O adapters:** parse trajectory/movie formats into a **canonical per-snapshot dictionary**.
- **Feature builders:** compute framework mesh, surfaces, maxima/windows/connections payloads from the canonical dictionary.
- **Scene model:** unified dictionary schema for traces/layout/camera/animation.
- **Renderer:** HTML/Three.js export from scene model.

Keep `devtools/scripts/md_pdb_threejs_viewer.py` as an example/CLI wrapper over the module, not the core implementation.

## 3) Unified Input Format (Canonical Snapshot Dictionary)

Each snapshot is represented as a dictionary so the same structure works for both fixed and flexible workflows:

```python
snapshot = {
    "atoms": ase.Atoms,           # framework + adsorbates for this frame
    "bonds": None | list,         # optional; if None -> auto-detect via update_with_covalent_bonds
    "cell": np.ndarray,           # 3x3 cell matrix (redundant if in atoms, but explicit for safety)
    "dgrid": None | np.ndarray,   # optional precomputed distance grid
    "regions": None | np.ndarray, # optional precomputed segmentation regions
    "connections": None | list,   # optional precomputed pore graph connections
    "meta": {...},                # arbitrary per-frame metadata (e.g., step number, energy)
}
```

- The I/O adapters (PDB, XYZ, DCD, etc.) produce a `list[snapshot]`.
- For **fixed framework**, only `atoms` (and optionally `bonds`) vary per frame; the first snapshot's `dgrid/regions/connections` are computed once and reused.
- For **flexible framework**, the user controls which frames trigger recomputation via `recompute_indices`.

## 4) Segmentation Computation Helper

A new helper in `mofography/movie/compute.py` (or `mofography/compute.py`) consumes the canonical snapshot dictionary and returns the segmentation payload:

```python
def compute_segmentation_for_snapshot(
    snapshot: dict,
    spacing: float = 0.25,
    mask_thickness: float = 0.5,
    h: float = 0.2,
    threshold_abs: float = 1.5,
    min_distance: float = 1.5,
    extra_padding: int = 5,
    use_rust: bool = True,
) -> dict:
    """
    Computes dgrid, regions, maxima, connections for a single snapshot.
    Returns a dict that can be merged back into the snapshot.
    """
```

- This keeps the segmentation logic in one place, reusable by both `visualize.py` and the movie module.
- It uses the existing `do_periodic_segmentation_rust` / `dgrid_from_atoms_cpu` pathways.

## 5) Lazy / Selective Computation Strategy

Because `dgrid` and `regions` are 3D arrays and can be large, **the default should NOT precompute them for every frame**.

- **Fixed framework (default):** compute segmentation once on the first snapshot, copy to all frames.
- **Flexible framework:** user provides `recompute_indices`:
  - `recompute_indices="all"` -> every frame (explicit opt-in, warned if >100 frames).
  - `recompute_indices=[0, 100, 200, ...]` -> only specified frames.
  - `recompute_indices=slice(None, None, 100)` -> every 100th frame.
- Uncomputed frames receive `None` for `dgrid/regions/connections`; the renderer skips rendering those traces for those frames (or falls back to the last computed frame).

## 6) Canonical Scene Model (Dictionary-First)

Top-level object:

```python
movie = {
    "version": "1.0",
    "frames": [...],
    "static": {...},
    "groups": {...},
    "styles": {...},
    "layout": {...},
    "meta": {...},
}
```

### 6.1 Frames

`frames[i]`:

```python
{
    "frame_index": int,
    "traces": [
        {
            "id": str,
            "type": "mesh" | "instanced_spheres" | "points" | "line_segments" | "unit_cell_lines",
            "name": str,
            "visible": bool,
            "style": {"color": ..., "opacity": ..., "roughness": ...},  # refs entry in top-level styles
            "data": {...},   # type-specific numeric payload
            "hover": [...] | None,  # per-vertex or per-instance hover strings (v1 pre-serialized; v2 templated)
        }
    ]
}
```

### 6.2 Static

Fixed traces that do not change per frame (e.g., framework mesh, unit cell, static pore graph).

### 6.3 Groups

Explicit group definitions so the JS renderer can auto-generate toggles without hardcoded linkage:

```python
"groups": {
    "internal_pore_graph": {
        "traces": ["internal_connections", "internal_nodes", "internal_maxima"],
        "label": "Internal Pore Graph"
    },
    "periodic_pore_graph": {
        "traces": ["periodic_connections", "periodic_nodes", "outside_maxima"],
        "label": "Periodic Pore Graph"
    },
    "guests": {
        "traces": ["adsorbates"],
        "label": "Adsorbates"
    }
}
```

### 6.4 Styles

Shared style registry (`DEFAULT_STYLES` in `mofography/movies.py`) to prevent drift between tooling and the Three.js export.

**Implemented behavior:** each trace carries a **resolved** `style` object (merged from the registry at build time). The HTML renderer uses that object directly (with string-key fallback for older payloads). Defaults are aligned with `devtools/scripts/md_pdb_threejs_viewer.py` (red connections, dashed periodic edges, yellow windows, pocket-colored internal maxima).

See **`docs/flexible_framework_movie_html.md`** for a practical flexible-framework + component workflow and sample RASPA paths.

### 6.5 Layout / Camera

```python
"layout": {
    "scene": {
        "camera": {"eye": [x, y, z], "center": [x, y, z], "up": [x, y, z]},
        "camera_keyframes": None | [...]  # Phase 2: optional per-frame camera interpolation
    }
}
```

## 7) API Proposal (Plotly-Like Ergonomics)

Core class:

```python
class MovieFigure:
    def __init__(self, data=None, layout=None, frames=None, groups=None, styles=None)
    def add_trace(trace, frame=None)
    def add_traces(traces, frame=None)
    def update_traces(patch, selector=None, frame=None)
    def update_layout(patch)
    def set_camera(eye=None, center=None, up=None, keyframes=None)
    def to_dict() -> dict
    @classmethod
    def from_dict(cls, d: dict) -> "MovieFigure"
    def to_html(path, template="threejs_default", sidecar=False)
```

Helper constructors (module-level functions):

```python
def make_framework_trace(atoms, bonds="auto", style="framework", name="Framework")
def make_surface_trace(dgrid, regions, atoms, isovalue=1.0, step_size=2, style="surface", name="Surface")
def make_maxima_trace(seg, atoms, style="internal_maxima", name="Maxima")
def make_windows_trace(seg, style="internal_windows", name="Windows")
def make_connections_trace(seg, atoms, style="internal_connections", name="Connections")
def make_guest_trace(atoms, style="guests", name="Adsorbates")
```

## 8) Reuse Plan for Existing Code

- **`mofography.analyze.update_with_covalent_bonds`**: the default `auto` bond mode. Replicates across PBC images as needed, stores bonds in `atoms.info['covalent_bonds']`.
- **`mofography.compute.do_periodic_segmentation_rust` / `dgrid_from_atoms_cpu`**: wrapped by the new `compute_segmentation_for_snapshot` helper.
- **`mofography.visualize` mesh builders** (`mesh_for_atoms`, etc.): keep for Plotly exports; the movie module uses its own numeric serializers (`_serialize_trimesh`, `_serialize_maxima_trace`, etc.) fed by the same raw numpy arrays.

**Decision**: Do not reuse Plotly trace constructors directly for the Three.js renderer. They are optimized for Plotly's JSON schema, not Three.js `BufferGeometry`. Instead, extract pure numpy computation into `compute_segmentation_for_snapshot`, then have separate thin serializers for Plotly vs. Three.js.

## 9) Supported Workflows

### A) Fixed Framework (Current Priority)

1. Parse trajectory -> `list[snapshot]`.
2. Compute bonds on first snapshot (or all if `bonds="auto"`).
3. Compute segmentation on first snapshot via `compute_segmentation_for_snapshot`.
4. Copy static traces (framework, surface, graph) to `movie["static"]`.
5. Build per-frame guest traces from `snapshot["atoms"]` positions.
6. Export HTML.

### B) Flexible Framework (Phase 2)

1. Parse trajectory -> `list[snapshot]`.
2. User specifies `recompute_indices` (explicit list, slice, or `"all"`).
3. Loop over `recompute_indices`, call `compute_segmentation_for_snapshot` for each.
4. Store computed segmentation back into the corresponding snapshots.
5. Build traces per frame; static traces are empty.
6. Export HTML (frames may be larger due to embedded meshes).

### C) Precompute + Render (Phase 3)

- Separate pipeline:
  - `preprocess_trajectory(trajectory_path, out_dir, recompute_indices=...)` -> writes cacheable artifacts (`.npz` for dgrid/regions, `.json` for trace metadata).
  - `render_from_artifacts(artifact_dir, out_html)` -> reads artifacts, builds scene model, emits HTML.
- Optional binary sidecar for large trajectories (see Section 11).

## 10) Input Format Expansion

- **v1 built-in adapters:**
  - RASPA PDB (current)
  - multi-frame XYZ / extXYZ (via `ase.io`)
  - ASE trajectory (`.traj`, `.db`)
- **v2 adapters (future):**
  - DCD / XTC / TRR (via `MDAnalysis` or `mdtraj`)
  - LAMMPS dump
- **Canonical per-frame ingest fields:**
  - unit cell (`a,b,c,alpha,beta,gamma` or full 3x3 matrix)
  - framework atom symbols + positions
  - adsorbate symbols + positions
  - optional bonds (list of pairs or `"auto"`)

## 11) Large Data / Streaming Strategy

Because a 10k-frame movie with embedded meshes can produce multi-hundred-MB HTML files, the renderer must support a sidecar mode:

- **Single-file mode (default):** all data embedded in one HTML as a JSON string. Suitable for < ~100 frames or small systems.
- **Sidecar mode (`sidecar=True` in `to_html`):**
  - Frames 0..N stored in a separate `.jsonl` or `.npz` file.
  - HTML contains only `static` traces + a frame index table.
  - JS renderer fetches frame chunks on demand via `fetch()`. This is essential for DCD/XTC trajectories that can be GBs.
- **Lazy loading in Python:** `MovieFigure` holds only frame metadata (indices, trace IDs) in memory; heavy arrays are loaded from disk when `to_html()` or `to_dict()` is called.

## 12) Renderer Template Mechanism

The `to_html(template="threejs_default")` method uses **Jinja2** templates stored in `mofography/movie/templates/`.

- `threejs_default.html`: the current r128-based viewer, parameterized by the scene model JSON.
- Templates reference `trace.data` arrays and `trace.style` keys from the scene model, eliminating hardcoded JS colors.
- `MANIFEST.in` must include `mofography/movie/templates/**/*`.

## 13) Hover Payload Serialization

- **v1 (current):** pre-serialize hover strings at build time (e.g., `MAXIMA.internal_hover`). Simple, but bloats JSON for many frames.
- **v2 optimization (Phase 3):** store structured hover dicts in the scene model and render client-side via a small JS template function. This reduces frame size by ~30-50% for text-heavy traces.

## 14) Milestones

- **M1:** finalize data model + `MovieFigure` skeleton + style registry
- **M2:** parity with current viewer traces (framework, surfaces, maxima/windows/connections, adsorbates) + fixed framework workflow
- **M3:** flexible-framework recompute support with `recompute_indices`
- **M4:** multi-format trajectory adapters (XYZ/extXYZ, ASE `.traj`)
- **M5:** precompute pipeline + sidecar renderer + optional parallel backend
- **M6:** docs/examples/tests and stabilization

## 15) Testing Strategy

- **Unit tests:**
  - trace schema validation (every trace must have `id`, `type`, `data`)
  - API mutation behavior (`add_trace`, `update_layout`, selectors)
  - adapter parsing correctness (PDB -> canonical snapshot)
  - `compute_segmentation_for_snapshot` parity with existing `do_periodic_segmentation_rust`
  - group-toggle consistency (internal/periodic linkage)
  - bond auto-detection fallback when `bonds=None`
- **Integration tests:**
  - fixed-framework HTML generation + JS smoke test (headless browser)
  - flexible-framework subset generation
  - hover payload integrity (maxima/windows/connections)
  - sidecar loading (fetch mock)
- **Regression tests:**
  - compare key counts/labels against trusted `visualize.py` outputs.
  - compare bond lists against `update_with_covalent_bonds` for a reference CIF.

## 16) Risks and Decisions

| Question | Decision |
|---|---|
| Reuse Plotly trace constructors directly? | **No.** Extract pure numpy computation, then have separate thin serializers for Plotly vs. Three.js. |
| Storage format for large trajectories? | **Sidecar `.jsonl`/`.npz` + JSON metadata.** Single-file default for small movies. |
| First-class trajectory formats in v1? | **RASPA PDB, multi-frame XYZ/extXYZ, ASE `.traj`.** DCD/XTC/TRR in v2. |
| Flexible framework: precompute upfront or lazy? | **Lazy.** Only compute for `recompute_indices`; uncomputed frames reuse last known or skip traces. |
| Compatibility between Plotly and Three.js defaults? | **Shared style registry (`styles` key in scene model).** Plotly and Three.js renderers both read from it. |
| Hover templating? | **Pre-serialized strings in v1; structured dicts + JS templates in v2.** |
