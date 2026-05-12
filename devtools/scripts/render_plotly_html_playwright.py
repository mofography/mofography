#!/usr/bin/env python3
"""
Render smooth frame-by-frame screenshots from a Plotly HTML using Playwright.

Workflow implemented from user spec:
- Start: framework + cell only
- Every 1.5 rotations: toggle/transition traces through workflow stages
- Camera rotates about z-axis with fixed eye-z, user-specified eye/up/center

This script captures PNG frames and can optionally assemble an MP4 with ffmpeg.
"""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
from pathlib import Path

from playwright.sync_api import sync_playwright


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Render Plotly HTML frames with Playwright.")
    p.add_argument("--html", type=Path, required=True, help="Path to Plotly HTML file.")
    p.add_argument("--out-dir", type=Path, default=Path("frames_playwright"), help="Output frames directory.")
    p.add_argument("--fps", type=float, default=30.0, help="Target video FPS (default: 30).")
    p.add_argument("--duration", type=float, default=30.0, help="Target duration in seconds (default: 30).")
    p.add_argument("--rotations-per-stage", type=float, default=1.5, help="Rotations between stage toggles.")
    p.add_argument(
        "--eye", nargs=3, type=float, metavar=("X", "Y", "Z"), default=(1.75, 0.0, 0.5),
        help="Camera eye (default: 1.75 0 0.5).",
    )
    p.add_argument(
        "--up", nargs=3, type=float, metavar=("X", "Y", "Z"), default=(0.0, 0.0, 1.0),
        help="Camera up (default: 0 0 1).",
    )
    p.add_argument(
        "--center", nargs=3, type=float, metavar=("X", "Y", "Z"), default=(0.0, 0.0, 0.0),
        help="Camera center (default: 0 0 0).",
    )
    p.add_argument("--width", type=int, default=1280, help="Viewport width (default: 1280).")
    p.add_argument("--height", type=int, default=720, help="Viewport height (default: 720).")
    p.add_argument("--channel", type=str, default="chrome", help="Playwright channel (default: chrome).")
    p.add_argument("--browser-path", type=str, default="", help="Optional explicit browser executable path.")
    p.add_argument("--mp4", type=Path, default=None, help="Optional output MP4 path to assemble with ffmpeg.")
    p.add_argument("--jpeg-quality", type=int, default=95, help="Reserved for future; PNG is used now.")
    return p


def _js_helpers() -> str:
    return r"""
() => {
  const gd = document.querySelector('.plotly-graph-div');
  if (!gd) throw new Error('No .plotly-graph-div found');

  function getTraces() {
    return gd.data.map((t, i) => ({i, name: (t.name || ''), type: (t.type || '')}));
  }

  async function setVisible(indices, value) {
    if (!indices || indices.length === 0) return;
    const vis = new Array(indices.length).fill(value);
    await Plotly.restyle(gd, {visible: vis}, indices);
  }

  async function setOpacity(indices, opacity) {
    if (!indices || indices.length === 0) return;
    const op = new Array(indices.length).fill(opacity);
    await Plotly.restyle(gd, {opacity: op}, indices);
  }

  async function setCamera(cam) {
    await Plotly.relayout(gd, {'scene.camera': cam});
  }

  return { ok: true, traces: getTraces() };
}
"""


def _match_indices(traces: list[dict], predicate) -> list[int]:
    out: list[int] = []
    for t in traces:
        try:
            if predicate(t):
                out.append(int(t["i"]))
        except Exception:
            continue
    return out


def _group_indices(traces: list[dict]) -> dict[str, list[int]]:
    def has(s: str):
        return lambda t: s.lower() in str(t.get("name", "")).lower()

    groups = {
        "framework": _match_indices(traces, lambda t: has("framework")(t)),
        "cell": _match_indices(traces, lambda t: has("unit cell")(t) or has("cell")(t)),
        "maxima_distinct": _match_indices(traces, lambda t: has("maxima (distinct)")(t)),
        "maxima_internal": _match_indices(traces, lambda t: has("internal maxima")(t)),
        "regions": _match_indices(traces, lambda t: has("regions")(t) and not has("pbc")(t)),
        "regions_pbc": _match_indices(traces, lambda t: has("regions pbc")(t) or has("pbc-merged")(t)),
        "regions_cluster": _match_indices(traces, lambda t: has("pocket types by geometry")(t) or has("clusters")(t)),
        "regions_cluster_trans": _match_indices(traces, lambda t: has("transparent")(t)),
        "internal_connections": _match_indices(traces, lambda t: has("internal connections")(t)),
        "windows": _match_indices(traces, lambda t: has("windows")(t)),
        "periodic_connections": _match_indices(traces, lambda t: has("periodic connections")(t)),
        "periodic_nodes": _match_indices(traces, lambda t: has("periodic maxima")(t)),
    }
    # De-duplicate cell by excluding framework accidental matches if needed
    if groups["framework"] and groups["cell"]:
        fw = set(groups["framework"])
        groups["cell"] = [i for i in groups["cell"] if i not in fw]
    return groups


def main() -> None:
    args = _build_parser().parse_args()
    html = args.html.resolve()
    if not html.is_file():
        raise SystemExit(f"HTML file not found: {html}")

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    fps = float(args.fps)
    duration = float(args.duration)
    total_frames = max(1, int(round(fps * duration)))
    rot_per_stage = float(args.rotations_per_stage)
    stage_thresholds = [rot_per_stage * i for i in range(1, 8)]  # 7 transitions
    total_rotations = stage_thresholds[-1]
    frames_per_rotation = total_frames / max(total_rotations, 1e-9)

    ex0, ey0, ez0 = map(float, args.eye)
    upx, upy, upz = map(float, args.up)
    cx, cy, cz = map(float, args.center)
    rxy = math.hypot(ex0, ey0)
    phi0 = math.atan2(ey0, ex0)

    print(
        f"Rendering {total_frames} frames ({duration:.1f}s @ {fps:.1f} fps), "
        f"~{frames_per_rotation:.1f} frames/rotation"
    )

    with sync_playwright() as p:
        launch_kwargs = dict(headless=True)
        if args.browser_path:
            launch_kwargs["executable_path"] = args.browser_path
        elif args.channel:
            launch_kwargs["channel"] = args.channel
        browser = p.chromium.launch(**launch_kwargs)
        page = browser.new_page(viewport={"width": args.width, "height": args.height})
        page.goto(html.as_uri(), wait_until="networkidle")
        page.wait_for_selector(".plotly-graph-div", timeout=120000)
        plot = page.locator(".plotly-graph-div")

        # Get trace metadata
        ret = page.evaluate(_js_helpers())
        traces = ret["traces"]
        groups = _group_indices(traces)
        all_idx = [int(t["i"]) for t in traces]

        async_eval = page.evaluate

        # Remove any built-in Plotly animation UI from source HTML
        # so this renderer acts on a plain figure view.
        async_eval(
            """() => {
                const gd = document.querySelector('.plotly-graph-div');
                return Plotly.relayout(gd, {
                    sliders: [],
                    updatemenus: []
                });
            }"""
        )

        # Initial visibility state: framework + cell only
        async_eval(
            """([allIdx, showIdx]) => {
                const gd = document.querySelector('.plotly-graph-div');
                const hide = allIdx.filter(i => !showIdx.includes(i));
                return Promise.all([
                    Plotly.restyle(gd, {visible: new Array(showIdx.length).fill(true)}, showIdx),
                    Plotly.restyle(gd, {visible: new Array(hide.length).fill('legendonly')}, hide),
                ]);
            }""",
            [all_idx, groups["framework"] + groups["cell"]],
        )
        # Ensure framework fully opaque at start
        async_eval(
            """(idx) => {
                const gd = document.querySelector('.plotly-graph-div');
                return Plotly.restyle(gd, {opacity: new Array(idx.length).fill(1.0)}, idx);
            }""",
            groups["framework"],
        )

        # Stage transition callbacks (cumulative visibility)
        # Requested order:
        # 0 start: framework + cell
        # 1 + local maxima + internal maxima
        # 2 + segmented surface (regions)
        # 3 + segmented surface (regions_pbc)
        # 4 + segmented surface (regions_cluster)
        # 5 + windows
        # 6 + internal connections
        # 7 + periodic connections + periodic nodes
        fired = [False] * len(stage_thresholds)

        def fire_transition(k: int) -> None:
            gd = ".plotly-graph-div"
            if k == 0:
                # Stage 2: local maxima + internal maxima (framework stays opaque)
                async_eval(
                    """([fw, toShow]) => {
                        const gd = document.querySelector('.plotly-graph-div');
                        return Promise.all([
                            Plotly.restyle(gd, {opacity: new Array(fw.length).fill(1.0)}, fw),
                            Plotly.restyle(gd, {visible: new Array(toShow.length).fill(true)}, toShow),
                        ]);
                    }""",
                    [groups["framework"], groups["maxima_distinct"] + groups["maxima_internal"]],
                )
            elif k == 1:
                # Stage 3: segmented surface (regions)
                async_eval(
                    """(idx) => {
                        const gd = document.querySelector('.plotly-graph-div');
                        return Plotly.restyle(gd, {visible: new Array(idx.length).fill(true)}, idx);
                    }""",
                    groups["regions"],
                )
            elif k == 2:
                # Stage 4: segmented surface (regions_pbc)
                async_eval(
                    """(idx) => {
                        const gd = document.querySelector('.plotly-graph-div');
                        return Plotly.restyle(gd, {visible: new Array(idx.length).fill(true)}, idx);
                    }""",
                    groups["regions_pbc"],
                )
            elif k == 3:
                # Stage 5: segmented surface (regions_cluster)
                async_eval(
                    """(idx) => {
                        const gd = document.querySelector('.plotly-graph-div');
                        return Plotly.restyle(gd, {visible: new Array(idx.length).fill(true)}, idx);
                    }""",
                    groups["regions_cluster"],
                )
            elif k == 4:
                # Stage 6: show windows, make clusters transparent for readability
                async_eval(
                    """([opaque, trans, windows]) => {
                        const gd = document.querySelector('.plotly-graph-div');
                        return Promise.all([
                            Plotly.restyle(gd, {visible: new Array(opaque.length).fill('legendonly')}, opaque),
                            Plotly.restyle(gd, {visible: new Array(trans.length).fill(true)}, trans),
                            Plotly.restyle(gd, {visible: new Array(windows.length).fill(true)}, windows),
                        ]);
                    }""",
                    [groups["regions_cluster"], groups["regions_cluster_trans"], groups["windows"]],
                )
            elif k == 5:
                # Stage 7: show internal connections
                async_eval(
                    """(idx) => {
                        const gd = document.querySelector('.plotly-graph-div');
                        return Plotly.restyle(gd, {visible: new Array(idx.length).fill(true)}, idx);
                    }""",
                    groups["internal_connections"],
                )
            elif k == 6:
                # Stage 8: show periodic connections and periodic maxima
                async_eval(
                    """([pc, pn]) => {
                        const gd = document.querySelector('.plotly-graph-div');
                        return Promise.all([
                            Plotly.restyle(gd, {visible: new Array(pc.length).fill(true)}, pc),
                            Plotly.restyle(gd, {visible: new Array(pn.length).fill(true)}, pn),
                        ]);
                    }""",
                    [groups["periodic_connections"], groups["periodic_nodes"]],
                )

        for i in range(total_frames):
            rotation = i / frames_per_rotation
            # stage transitions
            for k, thr in enumerate(stage_thresholds):
                if (not fired[k]) and rotation >= thr:
                    fire_transition(k)
                    fired[k] = True

            # camera update (z-axis rotation)
            phi = phi0 + (rotation * 2.0 * math.pi)
            cam = {
                "eye": {"x": math.cos(phi) * rxy, "y": math.sin(phi) * rxy, "z": ez0},
                "up": {"x": upx, "y": upy, "z": upz},
                "center": {"x": cx, "y": cy, "z": cz},
            }
            async_eval(
                """(cam) => {
                    const gd = document.querySelector('.plotly-graph-div');
                    return Plotly.relayout(gd, {'scene.camera': cam});
                }""",
                cam,
            )
            plot.screenshot(path=str(out_dir / f"frame_{i:04d}.png"))
            if i % 50 == 0 or i == total_frames - 1:
                print(f"[{i+1}/{total_frames}]")

        browser.close()

    if args.mp4:
        if shutil.which("ffmpeg") is None:
            raise SystemExit("ffmpeg not found on PATH; cannot assemble MP4.")
        mp4_path = args.mp4.resolve()
        mp4_path.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            "ffmpeg",
            "-y",
            "-framerate",
            str(fps),
            "-i",
            str(out_dir / "frame_%04d.png"),
            "-c:v",
            "libx264",
            "-pix_fmt",
            "yuv420p",
            str(mp4_path),
        ]
        subprocess.run(cmd, check=True)
        print(f"Wrote MP4: {mp4_path}")


if __name__ == "__main__":
    main()

