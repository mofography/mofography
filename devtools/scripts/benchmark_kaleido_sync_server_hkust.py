#!/usr/bin/env python3
"""
Benchmark Kaleido ``write_image`` x3 on the full HKUST workflow figure:

1) Without ``kaleido.start_sync_server()`` (default / per-call behaviour on v1).
2) After ``kaleido.start_sync_server()`` (persistent sync server; see Kaleido #400).

Requires the same env as ``generate_hkust_workflow_mp4.py`` (mofography, kaleido v1, …).
"""

from __future__ import annotations

import importlib.util
import shutil
import sys
import tempfile
import time
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent.parent


def _load_hkust_prepare():
    repo = _repo_root()
    path = repo / "generate_hkust_workflow_mp4.py"
    spec = importlib.util.spec_from_file_location("hkust_workflow_mp4_gen", path)
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(repo))
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def _time_three_writes(fig, w: int, h: int, tmp: Path, label: str) -> list[float]:
    times: list[float] = []
    for i in range(3):
        p = tmp / f"{label}_{i}.png"
        t0 = time.perf_counter()
        fig.write_image(str(p), width=w, height=h, scale=1)
        times.append(time.perf_counter() - t0)
    return times


def main() -> None:
    import importlib.metadata as im
    import plotly

    print("plotly:", plotly.__version__)
    try:
        import kaleido

        try:
            print("kaleido:", im.version("kaleido"))
        except Exception:
            print("kaleido:", getattr(kaleido, "__version__", kaleido))
    except ImportError as e:
        sys.exit(f"kaleido not importable: {e}")

    if not hasattr(kaleido, "start_sync_server"):
        sys.exit(
            "This benchmark needs kaleido v1+ with start_sync_server(). "
            "Upgrade kaleido / plotly per Plotly static export docs."
        )

    repo = _repo_root()
    hk = _load_hkust_prepare()

    print("Building HKUST workflow figure (same traces as MP4 export)…")
    t_build = time.perf_counter()
    ctx = hk.prepare_hkust_workflow_export(repo)
    print(f"  build wall time: {time.perf_counter() - t_build:.1f} s")

    fig = ctx["fig"]
    w, h = ctx["render_w"], ctx["render_h"]
    center = ctx["center_cart"]
    fig.update_layout(
        scene_camera=hk._camera_from_orbit(
            0.7, center, x_bias=ctx["x_bias"], orbit_r=ctx["orbit_r"]
        )
    )

    tmp = Path(tempfile.mkdtemp(prefix="kaleido_bench_"))
    try:
        kaleido.stop_sync_server(silence_warnings=True)

        print("\n--- A) Three write_image calls WITHOUT start_sync_server ---")
        t_a = _time_three_writes(fig, w, h, tmp, "noserver")

        print("\n--- B) start_sync_server(), then three write_image calls ---")
        kaleido.start_sync_server(silence_warnings=True)
        t_b = _time_three_writes(fig, w, h, tmp, "server")
        kaleido.stop_sync_server(silence_warnings=True)

        def fmt_ms(xs: list[float]) -> str:
            return ", ".join(f"{x * 1000:.0f} ms" for x in xs)

        print("\n=== Summary (per image) ===")
        print(f"  Without server: {fmt_ms(t_a)}  |  total {sum(t_a):.2f} s")
        print(f"  With server:    {fmt_ms(t_b)}  |  total {sum(t_b):.2f} s")
        ratio = sum(t_a) / max(sum(t_b), 1e-9)
        print(f"  Speedup (total B vs A): {ratio:.2f}x faster with server" if ratio > 1 else f"  Ratio A/B: {ratio:.2f}")
        print(f"\nPNG scratch dir: {tmp}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
