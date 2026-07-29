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

        page.evaluate("() => applyLineWidth('connections_internal', 2.5)")
        page.wait_for_timeout(200)
        pos_lw = conn_pos()
        assert pos_lw is not None
        dist2 = (pos_lw["x"] ** 2 + pos_lw["y"] ** 2 + pos_lw["z"] ** 2) ** 0.5
        assert abs(dist2 - dist0) < 1e-4, (
            f"pore graph jumped after line-width change: before={pos_exploded} after={pos_lw}"
        )

        browser.close()

    print("OK: exploded connections_internal position preserved after filter and line-width changes")
    return 0


def _assert_repr_visibility(page, repr_val: str, expected_visible: list[str], expected_hidden: list[str]) -> None:
    page.select_option("#reprSelect", repr_val)
    page.wait_for_timeout(150)
    state = page.evaluate(
        """(ids) => {
          const out = {};
          ids.forEach(id => {
            const o = objectsById[id];
            out[id] = !!(o && o.visible);
          });
          return out;
        }""",
        expected_visible + expected_hidden,
    )
    for tid in expected_visible:
        assert state.get(tid) is True, f"{tid} should be visible for {repr_val}, got {state}"
    for tid in expected_hidden:
        assert state.get(tid) is False, f"{tid} should be hidden for {repr_val}, got {state}"


def test_framework_repr_switching() -> int:
    path = HTML
    if not path.is_file():
        print(f"Missing {path}")
        return 1
    url = path.resolve().as_uri()
    ids = [
        "framework_atoms", "framework_bonds",
        "framework_stick_atoms", "framework_stick_bonds",
        "framework_spacefill_atoms",
    ]
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto(url, wait_until="networkidle", timeout=120_000)
        page.wait_for_function("typeof objectsById !== 'undefined'")

        _assert_repr_visibility(
            page,
            "ball_and_stick",
            ["framework_atoms", "framework_bonds"],
            ["framework_stick_atoms", "framework_stick_bonds", "framework_spacefill_atoms"],
        )
        _assert_repr_visibility(
            page,
            "stick",
            ["framework_stick_atoms", "framework_stick_bonds"],
            ["framework_atoms", "framework_bonds", "framework_spacefill_atoms"],
        )
        _assert_repr_visibility(
            page,
            "space_filling",
            ["framework_spacefill_atoms"],
            ["framework_atoms", "framework_bonds", "framework_stick_atoms", "framework_stick_bonds"],
        )

        page.evaluate("() => { explodeSlider.value = '1.5'; applyExplodeOffsets(); }")
        page.select_option("#reprSelect", "stick")
        page.wait_for_timeout(150)
        pos = page.evaluate(
            """() => {
              const o = objectsById.framework_stick_atoms;
              return o ? { visible: o.visible, x: o.position.x } : null;
            }"""
        )
        assert pos and pos["visible"] and abs(pos["x"]) > 1.0
        browser.close()

    print("OK: framework stick/spacefill repr switching and exploded offset")
    return 0


if __name__ == "__main__":
    code = main()
    if code != 0:
        raise SystemExit(code)
    raise SystemExit(test_framework_repr_switching())
