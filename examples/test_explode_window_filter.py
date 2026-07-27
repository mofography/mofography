#!/usr/bin/env python3
"""Playwright smoke test: exploded pore graph stays offset after min-window filter."""
from __future__ import annotations

import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
HTML = ROOT / "examples/viz_outputs/ZIF-8-RASPA_spacing0.25_threejs.html"


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else HTML
    if not path.is_file():
        print(f"Missing {path}")
        return 1

    url = path.resolve().as_uri()
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto(url, wait_until="networkidle", timeout=120_000)
        page.wait_for_function("typeof objectsById !== 'undefined' && objectsById.connections_internal")

        def conn_pos() -> dict:
            return page.evaluate(
                """() => {
                  const o = objectsById.connections_internal;
                  if (!o) return null;
                  return { x: o.position.x, y: o.position.y, z: o.position.z };
                }"""
            )

        page.evaluate("() => { explodeSlider.value = '1.5'; applyExplodeOffsets(); }")
        pos_exploded = conn_pos()
        assert pos_exploded is not None
        dist0 = (pos_exploded["x"] ** 2 + pos_exploded["y"] ** 2 + pos_exploded["z"] ** 2) ** 0.5
        assert dist0 > 1e-3, f"expected exploded offset, got {pos_exploded}"

        page.check("#hideEdgesBelowWindow")
        page.fill("#windowMinThreshold", "1.0")
        page.dispatch_event("#windowMinThreshold", "input")
        page.wait_for_timeout(200)

        pos_after = conn_pos()
        assert pos_after is not None
        dist1 = (pos_after["x"] ** 2 + pos_after["y"] ** 2 + pos_after["z"] ** 2) ** 0.5
        assert abs(dist1 - dist0) < 1e-4, (
            f"pore graph jumped after filter: before={pos_exploded} after={pos_after}"
        )

        periodic_visible = page.evaluate(
            "() => objectsById.connections_periodic && objectsById.connections_periodic.visible"
        )
        assert periodic_visible is False, "periodic connections should stay hidden"

        browser.close()

    print("OK: exploded connections_internal position preserved after min-window filter")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
