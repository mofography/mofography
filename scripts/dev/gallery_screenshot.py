#!/usr/bin/env python
"""Clean, canvas-only screenshots of a mofography Three.js scene (no control sidebar).

Drives the standalone HTML through its own DOM/JS (hide HUD, legend checkboxes,
Representation dropdown, Ortho, Exploded slider, view direction, Bg buttons) and
saves the WebGL canvas through the viewer's ``window.__mofFigureExportPng``.

Usage (needs internet: the page loads Three.js from a CDN):

    export PLAYWRIGHT_BROWSERS_PATH=$HOME/Library/Caches/ms-playwright
    PY=/Users/arungopalan/miniforge3/envs/mgr_env/bin/python

    # 1) see what legend entries the scene has
    $PY scripts/dev/gallery_screenshot.py docs/interactive/nu1000_pockets_windows_threejs.html --list

    # 2) take a picture
    $PY scripts/dev/gallery_screenshot.py docs/interactive/nu1000_pockets_windows_threejs.html \\
        docs/tutorial_images/gallery/nu1000_exploded_ortho.png \\
        --ortho --explode 0.6 --view 1 1 1 --only "Framework" "Pocket"

Legend options match by case-insensitive substring of the legend label.
  --only A B   switch every legend entry off except those matching A, B ...
  --on A       additionally switch matching entries on
  --off A      additionally switch matching entries off
  --zoom       camera distance multiplier, e.g. 1.3 to zoom out
  --repr       ball_and_stick | stick | space_filling
  --bg         "#ffffff" (solid) | "radial:#centre:#edge" | "linear:#c1:#c2" | transparent
  --jpeg Q     write JPEG (quality Q) instead of PNG (opaque backgrounds only)
"""
from __future__ import annotations

import argparse
import base64
import io
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

LEGEND_JS = """() => [...document.querySelectorAll('#groupChecks label')]
  .filter(l => l.querySelector('input[type=checkbox]'))
  .map(l => ({text: l.textContent.trim(), on: l.querySelector('input').checked}))"""

SET_LEGEND_JS = """([only, on, off]) => {
  const low = s => s.toLowerCase();
  const match = (t, keys) => keys.some(k => low(t).includes(low(k)));
  const labels = [...document.querySelectorAll('#groupChecks label')]
    .filter(l => l.querySelector('input[type=checkbox]'));
  const set = (cb, v) => { if (cb.checked !== v) cb.click(); };
  const missed = [];
  const isOpt = t => /screen-door|hide edges/i.test(t);  // options, not layers
  if (only.length) labels.forEach(l => { if (!isOpt(l.textContent)) set(l.querySelector('input'), match(l.textContent, only)); });
  labels.forEach(l => {
    const cb = l.querySelector('input');
    if (match(l.textContent, on)) set(cb, true);
    if (match(l.textContent, off)) set(cb, false);
  });
  [...only, ...on, ...off].forEach(k => {
    if (!labels.some(l => low(l.textContent).includes(low(k)))) missed.push(k);
  });
  return missed;
}"""

BG_JS = """(spec) => {
  const click = id => document.getElementById(id).click();
  if (spec === 'transparent') { click('bgGradientTransparent'); return; }
  const parts = spec.split(':');
  const kind = ['radial', 'linear'].includes(parts[0]) ? parts[0] : 'solid';
  const cols = kind === 'solid' ? [spec] : parts.slice(1);
  _bgSpec = {type: kind, colors: cols, active: true};
  document.getElementById('bgPicker').value = cols[0];
  renderer.setClearColor(0x000000, 1.0);
  setCheckerboardBackground(false);
  applyBgSpec();
}"""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("html")
    ap.add_argument("out", nargs="?")
    ap.add_argument("--list", action="store_true", help="print legend labels and exit")
    ap.add_argument("--width", type=int, default=1600)
    ap.add_argument("--height", type=int, default=1000)
    ap.add_argument("--only", nargs="*", default=[])
    ap.add_argument("--on", nargs="*", default=[])
    ap.add_argument("--off", nargs="*", default=[])
    ap.add_argument("--repr", dest="repr_", choices=["ball_and_stick", "stick", "space_filling"])
    ap.add_argument("--ortho", action="store_true")
    ap.add_argument("--explode", type=float, default=0.0)
    ap.add_argument("--view", nargs=3, type=float, metavar=("X", "Y", "Z"))
    ap.add_argument("--zoom", type=float, default=1.0,
                    help="camera distance multiplier (>1 zooms out; default views are cropped at 16:10)")
    ap.add_argument("--bg")
    ap.add_argument("--jpeg", type=int, metavar="Q")
    ap.add_argument("--wait", type=float, default=2.0, help="seconds to let the scene settle")
    ap.add_argument("--eval", dest="eval_", help="extra JS to run before the shot")
    a = ap.parse_args()
    if not a.list and not a.out:
        ap.error("out is required unless --list")

    url = Path(a.html).resolve().as_uri()
    with sync_playwright() as p:
        browser = p.chromium.launch(args=["--use-angle=swiftshader", "--enable-unsafe-swiftshader",
                                          "--ignore-gpu-blocklist"])
        page = browser.new_page(viewport={"width": a.width, "height": a.height})
        errs: list[str] = []
        page.on("pageerror", lambda e: errs.append(str(e)))
        page.goto(url, wait_until="load", timeout=120_000)
        page.wait_for_selector("#groupChecks label", timeout=120_000)
        page.wait_for_timeout(1500)

        if a.list:
            for e in page.evaluate(LEGEND_JS):
                print(("[x] " if e["on"] else "[ ] ") + e["text"])
            browser.close()
            return 0

        # hide the sidebar, then fire resize so renderer + camera use the full window
        page.evaluate("setHudVisible(false); window.dispatchEvent(new Event('resize'))")
        page.wait_for_timeout(300)
        missed = page.evaluate(SET_LEGEND_JS, [a.only, a.on, a.off])
        if missed:
            print("WARNING: no legend match for", missed, file=sys.stderr)
        if a.repr_:
            page.evaluate("v => { const s = document.getElementById('reprSelect');"
                          " s.value = v; s.dispatchEvent(new Event('change')); }", a.repr_)
        if a.ortho and not page.evaluate("camera.isOrthographicCamera"):
            page.evaluate("document.getElementById('toggleOrtho').click()")
        if a.explode:
            page.evaluate("v => { const s = document.getElementById('explodeSlider');"
                          " s.value = v; s.dispatchEvent(new Event('input')); }", a.explode)
        if a.view:
            page.evaluate("v => { ['X', 'Y', 'Z'].forEach((k, i) =>"
                          " document.getElementById('viewDir' + k).value = v[i]);"
                          " document.getElementById('applyViewDir').click(); }", a.view)
        if a.zoom != 1.0:
            page.evaluate("z => { camera.position.sub(controls.target).multiplyScalar(z).add(controls.target);"
                          " controls.update(); updateOrthoFrustum(); }", a.zoom)
        if a.bg:
            page.evaluate(BG_JS, a.bg)
        if a.eval_:
            page.evaluate(a.eval_)
        page.wait_for_timeout(int(a.wait * 1000))
        data = page.evaluate("window.__mofFigureExportPng(1)")
        browser.close()

    raw = base64.b64decode(data.split(",", 1)[1])
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if a.jpeg:
        from PIL import Image
        im = Image.open(io.BytesIO(raw)).convert("RGB")
        im.save(out, "JPEG", quality=a.jpeg, optimize=True)
    else:
        out.write_bytes(raw)
    print(f"wrote {out} ({out.stat().st_size/1024:.0f} KB)" + (f"  JS errors: {errs}" if errs else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
