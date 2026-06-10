# Orthographic Mode Issues — Handoff

## Context

We built a `Figure3D` API (Plotly-like wrapper around `MovieFigure`) in `mofography/movies.py` that renders static 3D scenes and MD movies via Three.js (r128). Two HTML files were generated for testing:

| File | Size | Content |
|------|------|---------|
| `cubtc_static.html` | 55 MB | CuBTC: framework + cell + 3 surfaces (raw/PBC/cluster) + maxima + connections + windows |
| `lta96_k_movie_v2.html` | 48 MB | LTA96: 50-frame K+ diffusion movie with framework + surface + pore graph + guests |

**Regenerate with:** `python scripts/dev/gen_cubtc_demo.py` (requires `mgr_env` conda environment).

**Run tests:** `pytest mofography/tests/test_movies.py -v` (30 tests, all passing).

---

## Issue 1: Orthographic lighting/coloring — surfaces look pale/washed out

### Current state

The HTML viewer has a "Toggle Orthographic" button in the HUD (next to pivot buttons). Clicking it switches between perspective (default) and orthographic cameras. What happens on toggle:

1. Camera switches (position and target preserved)
2. `updateOrthoFrustum()` computes frustum from camera-to-target distance: `frustumSize = 2 * dist * tan(22.5°)`. This makes the ortho frustum match the perspective visible area.
3. SSAO pass is **destroyed/recreated** — SSAO only works with perspective (produces black output with ortho). In ortho mode, `composer = null` and `renderer.render()` is used directly.
4. Directional lights get a **1.4× boost** to compensate for missing AO darkening.
5. Point materials get `depthTest: false` / `depthWrite: false` (see Issue 2).

### Problem

Without SSAO, all surfaces are evenly lit — no ambient occlusion depth shading. The 1.4× directional boost helps but doesn't fully compensate. The scene looks "pale" and "flat" compared to perspective mode.

### Possible approaches

**A. Baked per-vertex AO** (best quality, static scenes only)
- For static meshes (framework, surfaces), compute occlusion per vertex during mesh generation
- Ray-cast against neighboring atoms, count hits → darken vertex color
- Zero runtime cost, works with both projections
- ~50 lines of numpy in `gen_cubtc_demo.py`'s `make_surface_trace_smooth()` function

**B. Better SSAO alternative for ortho**
- Could try a different AO approach that works with ortho
- r128's SSAO pass is perspective-only; upgrading Three.js version might help
- Or implement a simple distance-based AO in a custom shader

**C. Material/lighting tweaks**
- Current directional boost: 1.4× (line ~2040 in movies.py toggle handler)
- Could try: higher boost, or adding a HemisphereLight intensity change
- Or adjust mesh material roughness dynamically for ortho

**D. Accept the difference**
- Ortho mode is for measurement/alignment, not for final renders
- Perspective mode (with SSAO) is for visual quality

### Where to look

- `mofography/movies.py` — template toggle handler at `toggleOrthoBtn.addEventListener` (~line 2020)
- `mofography/movies.py` — `recreateComposer()` function (~line 2000)
- `mofography/movies.py` — SSAO init (`if (!isOrtho)`) at ~line 2054
- `gen_cubtc_demo.py` — `make_surface_trace_smooth()` at top (~line 18)
- `mofography/visualize.py` — `default_mesh_lighting` dict for reference values
- `mofography/visualize.py` — `make_fig_orthographic()` for how Plotly handles it

---

## Issue 2: Scatter traces (maxima, windows) too small in orthographic mode

### Current state

Point markers (maxima, windows, connection nodes) use `THREE.PointsMaterial` with a circular sprite texture (`circleTex`, generated from a 64×64 canvas). Size multiplier: 1.3× (to compensate for circles appearing smaller than default squares). 

In ortho mode, points were invisible due to `depthTest` blocking them (linear depth in ortho puts them behind surfaces). Fixed: `depthTest: false, depthWrite: false` on toggle.

Current sizes:
- Maxima: 7.15 (world units)
- Windows: 3.25
- Connection nodes: 1.43

### Problem

User reports points appear smaller/correct in perspective but too small in orthographic. These sizes produce ~152px point diameter in ortho (at frustum 16.87, viewport 900px) — should be plenty visible. With `depthTest: false` they render on top of everything.

Possible remaining issues:
1. **No sizeAttenuation difference** — both modes use `sizeAttenuation: true`. In ortho, size is `value × (viewport / frustum)` = constant screen size regardless of distance. In perspective, size also depends on distance. The 1.3× circle compensation might not be enough in ortho.
2. **Color contrast** — without SSAO darkening, the bright scene makes colored points blend into surfaces. Maybe increase point opacity or border.
3. **The 1.3× multiplier is applied in `addPoints`** (`const size = ... * 1.3`) and `replacePointsGeometry` (`pointsBaseSize(trace) * 1.3`). Try increasing to 1.8× or 2.0× for ortho specifically.

### Where to look

- `mofography/movies.py` — `addPoints` function at ~line 2145, size line: `const size = (...) * 1.3;`
- `mofography/movies.py` — `replacePointsGeometry` at ~line 2248, size line: `pointsBaseSize(trace) * 1.3;`
- `mofography/movies.py` — circle texture creation at ~line 1955 (`ptCanvas`, `circleTex`)
- `mofography/movies.py` — toggle handler point material tweaks at ~line 2045

---

## Quick smoke test

```bash
# Regenerate both HTMLs
cd /Users/arungopalan/repos/mofography
/Users/arungopalan/miniforge3/envs/mgr_env/bin/python scripts/dev/gen_cubtc_demo.py

# Verify with Playwright
/Users/arungopalan/miniforge3/envs/mgr_env/bin/python -c "
from playwright.sync_api import sync_playwright; from pathlib import Path
with sync_playwright() as p:
    b=p.chromium.launch(headless=True); page=b.new_page(viewport={'width':1280,'height':900})
    page.goto(f'file://{Path(\"cubtc_static.html\").resolve()}',timeout=30000)
    page.wait_for_selector('canvas',timeout=15000);
    info=page.evaluate('''()=>{return{nObj:Object.keys(objectsById).length,ssao:!!composer,errs:[]};}''')
    page.click('#toggleOrtho');
    page.wait_for_timeout(2000);
    info2=page.evaluate('''()=>{return{ortho:!camera.isPerspectiveCamera,ssao:!!composer};}''')
    b.close()
    print(f'Persp: {info[\"nObj\"]}obj SSAO:{info[\"ssao\"]}  Ortho: SSAO:{info2[\"ssao\"]}')
"
```
