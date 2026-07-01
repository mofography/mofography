# Workflow Animation Handoff

This document is for handing off the Plotly workflow-animation task to another agent.

## Goal

Create a **single interactive Plotly HTML** (for browser playback + screen recording) that shows the pore-analysis workflow in ordered stages with the **legacy centered camera POV** (not orbital camera drift).

Primary target script:

- `devtools/scripts/workflow_animation_html.py`

Related scripts/context:

- `devtools/scripts/cif_camera_animation_html.py` (legacy centered z-spin camera proof-of-concept)
- `generate_hkust_workflow_mp4.py` (same analysis stack, MP4 path via Kaleido)
- `mofography/visualize.py` (all plotting functions/traces used below)

---

## Where To Look In Code

### Plot-building functions (in `mofography/visualize.py`)

- `plot_framework_as_mesh`
- `plot_cell`
- `plot_maxima`
- `plot_surface_dgrid`
- `plot_windows`
- `plot_connections_as_it_is`
- (optional alternatives) `get_isosurface_submeshes`, `plot_percolated_cluster_meshes`

### Analysis/segmentation helpers used in script

- `dgrid_from_atoms_rust`
- `do_periodic_segmentation_rust`
- `rag_from_connections`
- `add_maxima_to_rag`
- `add_pbc_to_rag`
- `add_psd_to_rag_with_pbc`
- `add_clustering_to_rag_geometry_js`
- `clustering_grid_from_regions`

---

## Camera / Motion Requirements (must keep)

Use the **original centered look**:

- Eye rotates in xy plane, fixed z
- `center = (0, 0, 0)`
- `up = (0, 0, 1)`
- No “planet orbit” around off-origin center

Reference camera function:

- in `devtools/scripts/workflow_animation_html.py`: `_camera_legacy_z_spin()`

Recommended defaults:

- `radius = 2.0`
- `z = 1.2`
- `turns_per_stage = 1.5`
- `n_frames_stage = 24`
- `frame_ms = 45`

---

## Surface Quality Requirement

Set distance-grid resolution (spacing) to:

- `spacing = 0.25` (default)

Reason: coarse spacing makes transparency/lighting look poor.

---

## Stage-by-Stage Trace Plan (requested behavior)

Use these workflow scenes in order:

1. **Framework + cell**
   - `plot_framework_as_mesh`
   - `plot_cell`

2. **Add maxima (before segmentation)**
   - Stage 1 traces +
   - `plot_maxima`

3. **Segmented surfaces shown opaque**
   - 3a: `plot_surface_dgrid(..., region_labels=regions, opacity=1.0)`
   - 3b: `plot_surface_dgrid(..., region_labels=regions_pbc-1, opacity=1.0)`
   - 3c: `plot_surface_dgrid(..., region_labels=regions_cluster-1, opacity=1.0)`

4. **Windows on transparent/dimmed surface**
   - transparent cluster surface (copy of 3c trace, lower opacity)
   - `plot_windows`

5. **Pore graph on transparent surface (internal)**
   - same transparent surface
   - `plot_connections_as_it_is(..., plot_periodic=False)`

6. **Periodic connections emphasized**
   - same transparent surface
   - `plot_connections_as_it_is(..., plot_periodic=True, emphasized line styling)`

7. **Final complete picture**
   - framework + cell + maxima
   - opaque cluster surface
   - windows
   - periodic graph

### Framework visibility note

When surfaces appear (Stages 3a/3b/3c), framework can be hidden for clarity.  
Final stage should include framework again with opaque surface.

---

## Current Script Status

`devtools/scripts/workflow_animation_html.py` currently implements:

- stage order above
- opaque -> transparent transition
- periodic-emphasis stage
- final full stage
- legacy centered z-spin camera
- slider + Play/Pause controls

---

## Run Commands

### Full-quality HTML (recommended for recording)

```bash
cd /Users/arungopalan/repos/mofography
python devtools/scripts/workflow_animation_html.py \
  mofography/examples/TestCIFs/P1-Cu-BTC.cif \
  -o /tmp/hkust_workflow_animation.html \
  --spacing 0.25 \
  --n-frames-stage 24 \
  --turns-stage 1.5 \
  --frame-ms 45 \
  --cdn
```

### Fast smoke test

```bash
python devtools/scripts/workflow_animation_html.py \
  mofography/examples/TestCIFs/P1-Cu-BTC.cif \
  -o /tmp/hkust_workflow_quick.html \
  --spacing 0.5 \
  --n-frames-stage 4 \
  --turns-stage 0.5 \
  --frame-ms 40 \
  --cdn
```

---

## Acceptance Criteria

- Camera feels like first good version: centered object, stable z-up spin.
- Stage transitions match requested sequence exactly.
- Segmentation surfaces are opaque when introduced, transparent for windows/graph.
- Final stage displays everything together clearly.
- `spacing=0.25` works without crashes and gives smoother surfaces.
- HTML loads and plays in browser for screen recording.

