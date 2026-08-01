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
                  const grp = objectsById.pore_explode_group;
                  if (grp) return { x: grp.position.x, y: grp.position.y, z: grp.position.z };
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
              const grp = objectsById.framework_explode_group;
              return o ? {
                visible: o.visible,
                grpX: grp ? grp.position.x : 0,
              } : null;
            }"""
        )
        assert pos and pos["visible"] and abs(pos["grpX"]) > 1.0
        browser.close()

    print("OK: framework stick/spacefill repr switching and exploded offset")
    return 0


def test_exploded_cell_boxes() -> int:
    path = HTML
    if not path.is_file():
        print(f"Missing {path}")
        return 1
    url = path.resolve().as_uri()
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto(url, wait_until="networkidle", timeout=120_000)
        page.wait_for_function("typeof objectsById !== 'undefined'")

        page.evaluate("() => { explodeSlider.value = '2.0'; applyExplodeOffsets(); }")
        state = page.evaluate(
            """() => {
              function cellState(id) {
                const o = objectsById[id];
                return o ? { visible: o.visible, pos: o.position.toArray() } : null;
              }
              const fwGrp = objectsById.framework_explode_group;
              const pbcGrp = objectsById.pbc_explode_group;
              const ar = objectsById.aromatic_rings;
              const fw = objectsById.framework_atoms;
              return {
                cell_framework: cellState('cell_framework'),
                cell_copy: cellState('cell_copy'),
                cell_pbc: cellState('cell_pbc'),
                fw_grp_children: fwGrp ? fwGrp.children.length : 0,
                pbc_grp_children: pbcGrp ? pbcGrp.children.length : 0,
                fw_parent: fw && fw.parent ? fw.parent.name : null,
                ar_parent: ar && ar.parent ? ar.parent.name : null,
                fw_pos: fw ? fw.position.toArray() : null,
                ar_pos: ar ? ar.position.toArray() : null,
              };
            }"""
        )
        assert state["cell_framework"] and state["cell_framework"]["visible"]
        assert abs(state["cell_framework"]["pos"][0]) > 5.0
        assert state["cell_copy"] and state["cell_copy"]["visible"]
        assert abs(state["cell_copy"]["pos"][0]) > 5.0
        assert state["cell_pbc"] and state["cell_pbc"]["visible"]
        assert state["cell_pbc"]["pos"][1] > 5.0
        assert state["fw_grp_children"] >= 3
        assert state["pbc_grp_children"] >= 1
        assert state["fw_parent"] == "framework_explode_group"
        assert state["ar_parent"] == "framework_explode_group"
        assert state["fw_pos"] == state["ar_pos"]
        browser.close()

    print("OK: exploded cell boxes visible and framework/aromatic grouped")
    return 0


def test_aromatic_rings_rendered() -> int:
    path = HTML
    if not path.is_file():
        print(f"Missing {path}")
        return 1
    url = path.resolve().as_uri()
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto(url, wait_until="networkidle", timeout=120_000)
        page.wait_for_function("typeof objectsById !== 'undefined'")

        state = page.evaluate(
            """() => {
              const rings = objectsById.aromatic_rings;
              const pent = objectsById.aromatic_highlights_pent;
              return {
                rings: rings ? { visible: rings.visible, count: rings.count } : null,
                pent: pent ? { visible: pent.visible, count: pent.count } : null,
              };
            }"""
        )
        assert state["rings"] and state["rings"]["visible"] and state["rings"]["count"] == 36
        assert state["pent"] and state["pent"]["visible"] and state["pent"]["count"] == 36

        page.evaluate("() => { explodeSlider.value = '1.5'; applyExplodeOffsets(); }")
        pos = page.evaluate(
            """() => {
              const o = objectsById.aromatic_rings;
              const grp = objectsById.framework_explode_group;
              return o ? {
                visible: o.visible,
                x: o.position.x,
                parent: o.parent ? o.parent.name : null,
                grpX: grp ? grp.position.x : null,
              } : null;
            }"""
        )
        assert pos and pos["visible"] and pos["parent"] == "framework_explode_group"
        assert abs(pos["grpX"]) > 1.0
        browser.close()

    print("OK: aromatic ring instanced meshes render in viewer")
    return 0


if __name__ == "__main__":
    code = main()
    if code != 0:
        raise SystemExit(code)
    code = test_framework_repr_switching()
    if code != 0:
        raise SystemExit(code)
    code = test_exploded_cell_boxes()
    if code != 0:
        raise SystemExit(code)
    raise SystemExit(test_aromatic_rings_rendered())
