"""Export Plotly figures to MP4 by rendering PNG frames (Kaleido) and stitching with ffmpeg.

Kaleido v1+ normally launches a renderer per ``write_image`` call, which is very slow for
many frames. These functions wrap rendering in ``kaleido.start_sync_server()`` by default
so one persistent Chromium session is reused (see Plotly/Kaleido `#400` and related issues).

Plotly does not write video directly. Typical workflow::

    # Build one Figure with all traces; control visibility per stage.
    fig = go.Figure(data=all_traces, layout=layout)
    stages = [
        ("Stage 1", [True, False, ...], 20),   # title, visibility per trace, repeat count
        ("Stage 2", [True, True, ...], 20),
    ]
    export_staged_figure_to_mp4(fig, "out.mp4", stages=stages, fps=12)
    # Or per-stage camera: ``frame_hook_per_stage=...`` (see docstring).

Requires ``python-kaleido`` and ``ffmpeg`` on PATH.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager, nullcontext
from pathlib import Path
from typing import TypeAlias

import plotly.graph_objects as go
from tqdm import tqdm

StageSpec: TypeAlias = tuple[Sequence[bool], int] | tuple[str | None, Sequence[bool], int]


@contextmanager
def kaleido_persistent_session(
    *,
    workers: int | None = None,
    silence_warnings: bool = True,
) -> Iterator[None]:
    """Start Kaleido's sync server for fast repeated static export, then stop it.

    With Kaleido 1.x, call ``start_sync_server()`` once before many ``fig.write_image``
    (or ``plotly.io.write_image``) calls so the same engine is reused instead of spawning
    Chromium for every frame.

    Parameters
    ----------
    workers
        If set, passed as ``n=...`` to ``start_sync_server`` (parallel renderer tabs when
        supported by your Kaleido build).
    silence_warnings
        Passed through to Kaleido start/stop helpers.
    """
    started = False
    try:
        import kaleido

        if hasattr(kaleido, "start_sync_server"):
            kwargs: dict = {"silence_warnings": silence_warnings}
            if workers is not None:
                kwargs["n"] = workers
            kaleido.start_sync_server(**kwargs)
            started = True
        yield
    finally:
        if started:
            try:
                import kaleido

                if hasattr(kaleido, "stop_sync_server"):
                    kaleido.stop_sync_server(silence_warnings=silence_warnings)
            except Exception:
                pass


def _kaleido_session_cm(
    *,
    use_sync_server: bool,
    workers: int | None,
    silence_warnings: bool,
):
    if use_sync_server:
        return kaleido_persistent_session(
            workers=workers, silence_warnings=silence_warnings
        )
    return nullcontext()


def set_traces_visible(fig: go.Figure, visibility: Sequence[bool]) -> None:
    """Set ``visible`` on each trace in order."""
    if len(visibility) != len(fig.data):
        raise ValueError(
            f"visibility length {len(visibility)} != number of traces {len(fig.data)}"
        )
    for trace, vis in zip(fig.data, visibility):
        trace.visible = vis


def _parse_stage(stage: StageSpec) -> tuple[str | None, Sequence[bool], int]:
    if len(stage) == 2:
        return None, stage[0], stage[1]
    if len(stage) == 3:
        return stage[0], stage[1], stage[2]
    raise TypeError(f"Each stage must be (visibility, n_frames) or (title, visibility, n_frames); got {stage!r}")


def export_staged_figure_to_mp4(
    fig: go.Figure,
    output_mp4: str | Path,
    *,
    stages: Sequence[StageSpec],
    fps: int = 12,
    width: int = 1280,
    height: int = 720,
    scale: int = 1,
    frame_hook: Callable[[go.Figure, int], None] | None = None,
    frame_hook_per_stage: Callable[[go.Figure, int, int, int], None] | None = None,
    title_font_size: int = 20,
    crf: int = 20,
    frames_dir: str | Path | None = None,
    cleanup_frames: bool = True,
    show_progress: bool = True,
    kaleido_sync_server: bool = True,
    kaleido_sync_workers: int | None = None,
    kaleido_silence_warnings: bool = True,
) -> Path:
    """Render staged visibility sequences to an H.264 MP4.

    Parameters
    ----------
    fig
        Plotly figure (mutates trace ``visible`` flags; restore manually if reusing ``fig``).
    output_mp4
        Output ``.mp4`` path.
    stages
        Each entry is either ``(visibility, n_frames)`` or ``(title, visibility, n_frames)``.
        ``visibility`` is a bool per trace, in ``fig.data`` order.
    fps
        Frames per second in the output video.
    width, height, scale
        Passed to ``fig.write_image`` (Kaleido).
    frame_hook
        Optional ``callback(fig, global_frame_index)`` after visibility (and title) are set,
        e.g. to orbit ``fig.update_layout(scene_camera=...)``.
    frame_hook_per_stage
        If set, called instead of ``frame_hook``: ``callback(fig, stage_index, frame_in_stage, n_frames_in_stage)``
        for independent camera motion each stage (e.g. one or more full rotations per stage).
    title_font_size
        Font size for the stage title set from the stage tuple.
    crf
        libx264 quality (lower is higher quality; 18–23 is typical).
    frames_dir
        If set, PNGs are written here and ``cleanup_frames`` is treated as False for that directory.
    cleanup_frames
        If True and ``frames_dir`` is None, delete the temporary directory after ffmpeg succeeds.
    show_progress
        If True, show a tqdm bar while each Kaleido ``write_image`` frame is rendered.
    kaleido_sync_server
        If True (default), wrap PNG rendering in :func:`kaleido_persistent_session` so Kaleido
        v1+ reuses one Chromium-backed sync server across all ``write_image`` calls.
    kaleido_sync_workers
        Optional ``n=`` for ``kaleido.start_sync_server`` when your Kaleido build supports
        multiple parallel renderer tabs.
    kaleido_silence_warnings
        Passed to Kaleido start/stop helpers.
    """
    if shutil.which("ffmpeg") is None:
        raise RuntimeError("ffmpeg not found on PATH; install ffmpeg to assemble MP4.")

    output_mp4 = Path(output_mp4)
    output_mp4.parent.mkdir(parents=True, exist_ok=True)

    own_tmp = frames_dir is None
    if frames_dir is None:
        tmp = tempfile.mkdtemp(prefix="plotly_mp4_frames_")
        frames_path = Path(tmp)
    else:
        frames_path = Path(frames_dir)
        frames_path.mkdir(parents=True, exist_ok=True)

    if frame_hook is not None and frame_hook_per_stage is not None:
        raise ValueError("Pass only one of frame_hook and frame_hook_per_stage.")

    total_png = sum(_parse_stage(s)[2] for s in stages)
    frame_counter = 0
    try:
        with _kaleido_session_cm(
            use_sync_server=kaleido_sync_server,
            workers=kaleido_sync_workers,
            silence_warnings=kaleido_silence_warnings,
        ):
            with tqdm(
                total=total_png,
                desc="Kaleido PNG frames",
                unit="fr",
                disable=not show_progress,
            ) as pbar:
                for stage_index, stage in enumerate(stages):
                    title, visibility, n_frames = _parse_stage(stage)
                    set_traces_visible(fig, visibility)
                    if title is not None:
                        fig.update_layout(
                            title=dict(
                                text=title,
                                x=0.5,
                                xanchor="center",
                                font=dict(size=title_font_size, color="black"),
                            )
                        )
                    for frame_in_stage in range(n_frames):
                        if frame_hook_per_stage is not None:
                            frame_hook_per_stage(
                                fig, stage_index, frame_in_stage, n_frames
                            )
                        elif frame_hook is not None:
                            frame_hook(fig, frame_counter)
                        out_png = frames_path / f"frame_{frame_counter:04d}.png"
                        fig.write_image(
                            str(out_png), width=width, height=height, scale=scale
                        )
                        frame_counter += 1
                        pbar.update(1)

        if frame_counter == 0:
            raise ValueError("No frames rendered: `stages` produced zero repeats.")

        ffmpeg_cmd = [
            "ffmpeg",
            "-y",
            "-framerate",
            str(fps),
            "-i",
            str(frames_path / "frame_%04d.png"),
            "-c:v",
            "libx264",
            "-pix_fmt",
            "yuv420p",
            "-crf",
            str(crf),
            str(output_mp4),
        ]
        subprocess.run(ffmpeg_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        msg = e.stderr or e.stdout or str(e)
        raise RuntimeError(f"ffmpeg failed: {msg}") from e
    finally:
        if own_tmp and cleanup_frames and frames_path.exists():
            shutil.rmtree(frames_path, ignore_errors=True)

    return output_mp4


def export_plotly_frames_to_mp4(
    fig: go.Figure,
    output_mp4: str | Path,
    *,
    fps: int = 12,
    width: int = 1280,
    height: int = 720,
    scale: int = 1,
    include_initial: bool = True,
    repeats_per_frame: int = 1,
    crf: int = 20,
    frames_dir: str | Path | None = None,
    cleanup_frames: bool = True,
    frame_hook: Callable[[go.Figure, int], None] | None = None,
    kaleido_sync_server: bool = True,
    kaleido_sync_workers: int | None = None,
    kaleido_silence_warnings: bool = True,
) -> Path:
    """Export a figure that uses Plotly's ``frames`` API (slider / animation) to MP4.

    Frame patches are applied **cumulatively** in order, matching in-browser animation behavior.

    Parameters
    ----------
    include_initial
        If True, the first PNG is the figure's starting state before any ``fig.frames`` patch.
    repeats_per_frame
        Duplicate each rendered state this many times (slows motion between keyframes).
    kaleido_sync_server, kaleido_sync_workers, kaleido_silence_warnings
        Same meaning as in :func:`export_staged_figure_to_mp4`.
    """
    if not fig.frames:
        raise ValueError("Figure has no `frames`; use `export_staged_figure_to_mp4` instead.")

    if shutil.which("ffmpeg") is None:
        raise RuntimeError("ffmpeg not found on PATH; install ffmpeg to assemble MP4.")

    output_mp4 = Path(output_mp4)
    output_mp4.parent.mkdir(parents=True, exist_ok=True)

    own_tmp = frames_dir is None
    if frames_dir is None:
        tmp = tempfile.mkdtemp(prefix="plotly_mp4_frames_")
        frames_path = Path(tmp)
    else:
        frames_path = Path(frames_dir)
        frames_path.mkdir(parents=True, exist_ok=True)

    working = go.Figure(fig)
    working.frames = ()

    counter = [0]

    def write_shot() -> None:
        if frame_hook is not None:
            frame_hook(working, counter[0])
        out_png = frames_path / f"frame_{counter[0]:04d}.png"
        working.write_image(str(out_png), width=width, height=height, scale=scale)
        counter[0] += 1

    try:
        with _kaleido_session_cm(
            use_sync_server=kaleido_sync_server,
            workers=kaleido_sync_workers,
            silence_warnings=kaleido_silence_warnings,
        ):
            if include_initial:
                for _ in range(repeats_per_frame):
                    write_shot()

            for fr in fig.frames:
                for i, tr_update in enumerate(fr.data):
                    if i < len(working.data):
                        working.data[i].update(tr_update)
                if fr.layout is not None:
                    working.update_layout(fr.layout)
                for _ in range(repeats_per_frame):
                    write_shot()

        if counter[0] == 0:
            raise ValueError("No frames written.")

        ffmpeg_cmd = [
            "ffmpeg",
            "-y",
            "-framerate",
            str(fps),
            "-i",
            str(frames_path / "frame_%04d.png"),
            "-c:v",
            "libx264",
            "-pix_fmt",
            "yuv420p",
            "-crf",
            str(crf),
            str(output_mp4),
        ]
        subprocess.run(ffmpeg_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        msg = e.stderr or e.stdout or str(e)
        raise RuntimeError(f"ffmpeg failed: {msg}") from e
    finally:
        if own_tmp and cleanup_frames and frames_path.exists():
            shutil.rmtree(frames_path, ignore_errors=True)

    return output_mp4
