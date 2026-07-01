# Orthographic Mode Color Washout — Handoff

## Problem
When switching to orthographic projection in the Three.js viewer, surface colors look **pale / washed out** compared to perspective mode. The underlying vertex colors (Pastel palette) are correct — it's a rendering/lighting issue specific to orthographic cameras.

## Screenshot
- `~/Downloads/dull_ortho.png` — ortho mode showing pale/washed surfaces

## Root Cause (Likely)
In perspective mode, SSAO (`THREE.SSAOPass`) runs as a post-processing pass and adds ambient occlusion shadows in crevices, making colors look rich and 3D. In orthographic mode, SSAO produces a black screen (known Three.js r128 bug), so we skip it entirely (`composer = null`). Without SSAO, the flat uniform lighting makes vertex colors appear washed out.

## What We've Tried
| Approach | Result |
|---|---|
| Boost directional light 1.3×, reduce ambient to 0.07 | Slightly better but still pale |
| Stronger key (0.85), weaker fill (0.18), hemi 0.35 | Made it worse |
| Reduce roughness ×0.3, clearcoat 1.0, clearcoatRoughness 0.03, emissive 0.25 | Washed out colors further (too much white reflection) |
| Minimal 1.15× directional boost only, no material changes | Still pale |

## Key Files
- **`mofography/movies.py`** — main file containing the Three.js template
  - Template starts at `_DEFAULT_THREEJS_TEMPLATE` (~line 2250)
  - Ortho init block: search for `if (isOrtho)` after light setup (~line 2493)
  - Ortho toggle handler: search for `toggleOrthoBtn.addEventListener` (~line 2425)
  - SSAO setup: search for `SSAOPass` (~line 2500)
  - Material creation for meshes: search for `MeshPhysicalMaterial` in `addMesh` function (~line 2570)
  - Light definitions: search for `AmbientLight`, `HemisphereLight`, `DirectionalLight` (~line 2485)
- **`gen_cubtc_demo.py`** — generates `cubtc_static.html` test file
- **`cubtc_static.html`** — test output, open in browser

## Current Light Values (perspective baseline)
```js
ambient = new THREE.AmbientLight(0xffffff, 0.28);
hemiLight = new THREE.HemisphereLight(0xddeeff, 0x3b3b3b, 0.25);
dirKey = new THREE.DirectionalLight(0xffffff, 0.48);  // position (5, 10, 7.5)
dirFill = new THREE.DirectionalLight(0xffffff, 0.32); // position (-5, -9, -7.5)
```

## Current Ortho Overrides
```js
// Init and toggle: only 1.15× directional boost
dirKey.intensity *= 1.15;
dirFill.intensity *= 1.15;
// No material changes (roughness, clearcoat, emissive untouched)
```

## Current Material Defaults (from DEFAULT_STYLES)
```python
"surface": {
    "roughness": 0.15,
    "clearcoat": 0.8,
    "clearcoatRoughness": 0.05,
    "metalness": 0.05,
    "opacity": 1.0,
    "transparent": false,
}
```

## SSAO Config (perspective only)
```js
ssaoPass.kernelRadius = 0.8;
ssaoPass.kernelSize = 24;
ssaoPass.minDistance = 0.001;
ssaoPass.maxDistance = 0.05;
```

## Technical Constraints
- **Three.js r128** (CDN) — cannot upgrade without breaking other things
- SSAO `produces black in ortho` is a known r128 bug, cannot be fixed without upgrading
- Template uses `{{placeholder}}` substitution (not Jinja/f-strings)
- `DATA` JS variable holds entire scene dict (can be 50+ MB)
- MeshPhysicalMaterial used for Plotly-like PBR appearance
- Vertex colors are baked into geometry attribute `color` via `Float32BufferAttribute`

## Possible Approaches to Try
1. **Custom fake AO shader for ortho** — implement a screen-space darkening in corners/edges using `onBeforeCompile` on MeshPhysicalMaterial (but previous shader injection attempts broke rendering — see critical context below)
2. **Hemisphere light with stronger sky/ground contrast** — sky color bright, ground dark, to simulate directional shadow
3. **Second directional light from below** as a dim fill with warm color to add depth
4. **Post-processing fake AO for ortho** — use a simple darken pass that's not the broken SSAOPass
5. **Tone mapping** — apply `renderer.toneMapping = THREE.ACESFilmicToneMapping` with exposure adjustment for ortho
6. **Vertex color saturation boost** — in ortho mode, multiply vertex colors by a saturation factor > 1.0

## Critical Context — Shader Injection
Previous attempts to use `material.onBeforeCompile` to inject custom shader code **broke vertex color rendering** (patches of missing color on meshes). The issue was that `onBeforeCompile` in Three.js r128 can interfere with the vertex color pipeline in MeshPhysicalMaterial. **Any shader approach must be tested carefully** — generate the HTML and visually verify ALL surface regions show full color coverage with no missing patches.

## How to Test
```bash
cd /Users/arungopalan/repos/mofography
python scripts/dev/gen_cubtc_demo.py  # generates cubtc_static.html
open cubtc_static.html     # opens in browser
# Click "Ortho" button to switch, compare with perspective
```

## Run Tests
```bash
python -m pytest mofography/tests/test_movies.py -q
```
