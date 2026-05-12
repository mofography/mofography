# Figure3D Implementation Plan

## Goal

Add a Plotly-like Figure API to mofography that renders static 3D scenes via Three.js (fast WebGL, SSAO, sRGB, better performance than Plotly). Reuses the existing `MovieFigure` infrastructure and template from `mofography/movies.py`.

## Architecture (3-layer)

```
Layer 1: Shared compute functions (pure numpy → dicts)
  mofography/movies.py (promote existing _serialize_* to public)

Layer 2: Convenience builders on MovieFigure (chainable .add_* methods)

Layer 3: Jupyter integration (_repr_html_, show(), renderer='browser')
```

## Jupyter Integration

| Method | Implementation |
|--------|---------------|
| `_repr_html_()` | IFrame `srcdoc` pattern (pydeck proven). Avoids CSP issues across Jupyter environments. |
| `show(renderer='notebook')` | Calls `IPython.display.display(HTML(self._repr_html_()))` |
| `show(renderer='browser')` | Writes temp HTML → `webbrowser.open()` (like Plotly's `fig.show(renderer='browser')`) |

**Caveats:**
- Large scenes (50MB+) may be slow to embed in Jupyter — use `renderer='browser'` for those
- CDN scripts (three.js from jsdelivr) need internet access
- 2D charts (PSDs, histograms) stay Plotly — no Three.js equivalent needed

## API Design

```python
import mofography as mgr

# Build a static figure
fig = mgr.Figure3D()
fig.set_camera(eye=(2, 0, 1.35), center=(0, 0, 0), up=(0, 0, 1))

# Add traces (all chainable)
fig.add_framework(atoms, bonds='auto')
fig.add_cell(atoms)
fig.add_surface(dgrid, regions, atoms, isovalue=1.0, step_size=1)
fig.add_maxima(seg, atoms, hide_periodic=True)
fig.add_connections(seg, atoms, hide_periodic=True)
fig.add_windows(seg, hide_periodic=True)
fig.add_guests(atoms, color=[1.0, 0.41, 0.71])

# Export
fig.to_html("scene.html")      # file
fig.show()                      # Jupyter inline
fig.show(renderer='browser')    # open in browser
```

## What to build

### Phase 1: Promote compute functions (rename _serialize_* → compute_*)

File: `mofography/movies.py`

| Old name (private) | New name (public) | Returns |
|-------------------|-------------------|---------|
| `_serialize_framework_mesh` | `compute_framework_mesh` | `{positions, indices, colors, hover, n_vertices, n_triangles, n_atoms}` |
| `_serialize_cell_lines` | `compute_cell_lines` | `{positions, n_segments}` |
| `_serialize_segmentation_surface_mesh` | `compute_surface_mesh` | `{positions, indices, colors, labels, n_vertices, n_triangles}` |
| `_serialize_maxima_trace` | `compute_maxima_data` | `{internal_positions, internal_sizes, internal_colors, internal_hover, outside_positions, ...}` |
| `_serialize_windows_trace` | `compute_windows_data` | `{internal_positions, internal_hover, periodic_positions, ...}` |
| `_serialize_connection_traces` | `compute_connection_data` | `{internal_segments, internal_hover, periodic_segments, internal_points, ...}` |

Old names remain as aliases (backward compat). New API functions call the public versions.

### Phase 2: Add convenience methods to MovieFigure

File: `mofography/movies.py`

```python
class MovieFigure:
    # ... existing code ...

    def add_framework(self, atoms, bonds='auto', **style) -> 'MovieFigure'
    def add_cell(self, atoms, **style) -> 'MovieFigure'
    def add_surface(self, dgrid, regions, atoms, isovalue=1.0, step_size=1, **style) -> 'MovieFigure'
    def add_maxima(self, seg, atoms, hide_periodic=True, **style) -> 'MovieFigure'
    def add_connections(self, seg, atoms, hide_periodic=True, **style) -> 'MovieFigure'
    def add_windows(self, seg, hide_periodic=True, **style) -> 'MovieFigure'
    def add_guests(self, atoms, color=None, **style) -> 'MovieFigure'
```

Each method:
1. Calls the compute function
2. Wraps in `_make_trace`
3. Calls `self.add_trace()`
4. Returns `self` (chainable)

### Phase 3: Jupyter integration

File: `mofography/movies.py`

```python
class MovieFigure:
    def _repr_html_(self) -> str:
        """Jupyter display via IFrame srcdoc (avoids CSP issues)."""
        import html
        html_str = self._render_html_string()
        escaped = html.escape(html_str)
        return f'<iframe width="100%" height="600" frameborder="0" srcdoc="{escaped}"></iframe>'

    def _render_html_string(self) -> str:
        """Render to HTML string (without writing to file)."""
        ...

    def show(self, renderer='notebook') -> None:
        if renderer == 'browser':
            import tempfile, webbrowser
            path = Path(tempfile.mktemp(suffix='.html'))
            self.to_html(path)
            webbrowser.open(str(path))
        else:
            from IPython.display import display, HTML
            display(HTML(self._repr_html_()))
```

### Phase 4: High-level one-shot builder

File: `mofography/movies.py`

```python
def build_static_scene(
    atoms: Atoms,
    dgrid: np.ndarray = None,
    regions: np.ndarray = None,
    seg: dict = None,           # from do_periodic_segmentation_rust
    isovalue: float = 1.0,
    hide_periodic: bool = True,
    spacing: float = 0.25,
    **kwargs
) -> MovieFigure:
```

### Phase 5: Alias

```python
Figure3D = MovieFigure  # user-friendly alias
```

## Test Plan

### Existing test infrastructure

- `mofography/tests/test_movies.py` — already has MovieFigure tests (trace add/update, camera, template, guest spheres, HTML export). Add Figure3D-specific tests here.

### Regression tests (reference scripts to verify against)

| Script | CIF | Traces exercised | Weight |
|--------|-----|-----------------|--------|
| `test_cubtc.py` | P1-Cu-BTC.cif | surface, cell, framework, maxima, connections | **Primary** |
| `docs/generate_pore_graph_images.py` | P1-IRMOF-1.cif | surface, cell, framework, maxima, connections | **Primary** |
| `docs/test_percolated_cluster_coloring.py` | NU-1000_S1.cif | percolated cluster meshes | **Unique** |

### New test structure

Add to `mofography/tests/test_movies.py` (or new `test_figure3d.py`):

1. **test_smoke**: Load P1-IRMOF-1.cif, compute dgrid+segmentation, build Figure3D with all traces, export HTML, verify file exists + size > 0 + no JS errors
2. **test_jupyter_html**: Verify `_repr_html_()` returns valid IFrame wrapping `<html>` content
3. **test_browser_show**: `show(renderer='browser')` writes temp file, verify file content
4. **test_chainable**: Verify `fig.add_framework(...).add_cell(...).add_surface(...)` returns same fig instance
5. **test_non_regression**: Compare trace counts (n_vertices, n_triangles, n_segments) against golden values from test_cubtc.py output
6. **test_periodic_hidden**: Verify periodic traces have visible=False when hide_periodic=True

### Repo CIFs available for testing

```
mofography/examples/TestCIFs/P1-Cu-BTC.cif      (BTC)    ← primary test
mofography/examples/TestCIFs/P1-IRMOF-1.cif      (IRMOF)  ← primary test
mofography/examples/TestCIFs/NU-1000_S1.cif      (NU-1000)← percolation test
mofography/examples/TestCIFs/LTA96.cif           (LTA)    ← movie test
```

## Changes summary

| File | Changes | Lines |
|------|---------|-------|
| `mofography/movies.py` | Promote _serialize_* → compute_* (keep aliases) | ~30 |
| `mofography/movies.py` | 7 x .add_* convenience methods on MovieFigure | ~80 |
| `mofography/movies.py` | _repr_html_, _render_html_string, show() | ~50 |
| `mofography/movies.py` | build_static_scene() one-shot builder | ~50 |
| `mofography/movies.py` | Figure3D alias | 1 |
| `mofography/tests/test_movies.py` | 6 new test functions | ~150 |
| **Total** | | **~360** |

## Not in scope (v1)

- 2D chart migration (stays Plotly)
- Volume rendering (stays Plotly)
- ipywidgets backend (bidirectional Python↔JS communication)
- Camera keyframe animation (future feature)
- Modifying `visualize.py` (zero changes in v1)
