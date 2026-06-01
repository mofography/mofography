"""Plotly-parity helpers for granular Three.js scenes (Figure3D)."""
from __future__ import annotations
from typing import Any
import numpy as np
from ase import Atoms
from .movies import (
    Figure3D,
    MovieFigure,
    _make_trace,
    _serialize_connection_traces,
    make_cell_trace,
    make_connections_trace,
    make_framework_scatter_trace,
    make_maxima_trace,
    make_surface_trace,
    make_windows_trace,
)
__all__ = [
    "Figure3D", "MovieFigure", "Figure3JS", "qualitative_colorscale", "prepare_seg_threejs",
    "seg_from_connections", "add_trace", "add_traces", "plot_surface_dgrid",
    "plot_cell", "draw_framework", "plot_maxima", "plot_maxima_from_grid",
    "plot_connections", "plot_connections_as_it_is", "plot_windows",
]

def qualitative_colorscale(name: str = "Pastel") -> list[list[float]]:
    import plotly.express as px
    raw = getattr(px.colors.qualitative, name, None)
    if raw is None:
        raw = getattr(px.colors.sequential, name, px.colors.qualitative.Pastel)
    out: list[list[float]] = []
    for c in raw:
        s = str(c).strip()
        if s.startswith("rgb"):
            vals = [float(x) for x in s.replace("rgb(", "").replace("rgba(", "").replace(")", "").split(",")[:3]]
            out.append([v / 255.0 for v in vals])
        elif s.startswith("#"):
            h = s.lstrip("#")
            if len(h) == 3:
                h = "".join(ch * 2 for ch in h)
            out.append([int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4)])
        else:
            out.append([0.7, 0.7, 0.7])
    return out

def prepare_seg_threejs(seg: dict[str, Any], *, dgrid: np.ndarray | None = None) -> dict[str, Any]:
    out = dict(seg)
    if dgrid is not None:
        out["dgrid"] = np.asarray(dgrid)
    elif "dgrid" not in out:
        raise ValueError("seg must include dgrid or pass dgrid=")
    for key in ("regions", "connections", "all_maxima"):
        if key not in out:
            raise ValueError(f"seg missing required key {key!r}")
    return out

def seg_from_connections(connections, maxima, grid_shape, *, dgrid: np.ndarray | None = None) -> dict[str, Any]:
    shape = tuple(int(x) for x in grid_shape)
    if dgrid is None:
        dgrid = np.zeros(shape, dtype=np.float32)
    return {"dgrid": np.asarray(dgrid), "regions": np.zeros(shape, dtype=np.int32),
            "connections": list(connections), "all_maxima": np.asarray(maxima, dtype=float)}

def _hex_to_rgb01(color) -> list[float]:
    if isinstance(color, (list, tuple, np.ndarray)):
        arr = np.asarray(color, dtype=float).reshape(-1)
        if arr.size >= 3 and arr.max() <= 1.0:
            return [float(arr[0]), float(arr[1]), float(arr[2])]
        if arr.size >= 3:
            return [float(arr[0])/255, float(arr[1])/255, float(arr[2])/255]
    s = str(color).strip()
    if s.startswith("#"):
        h = s.lstrip("#")
        if len(h) == 3:
            h = "".join(ch*2 for ch in h)
        return [int(h[i:i+2], 16)/255 for i in (0,2,4)]
    return {"orange":[1,.65,0],"purple":[.5,0,.5],"yellow":[1,1,0],"red":[1,0,0]}.get(s.lower(), [.7,.7,.7])

def _style_from_plotly_kwargs(*, opacity=1.0, color=None, stylize=False, name="trace", **extra):
    style = {"opacity": float(opacity), "name": name}
    if color is not None:
        style["color"] = _hex_to_rgb01(color)
    if stylize:
        style.setdefault("metalness", 0.05)
        style.setdefault("roughness", 0.45)
    style.update(extra)
    return style

def add_trace(fig: Figure3D, trace: dict[str, Any]) -> Figure3D:
    fig.add_trace(trace)
    return fig

def add_traces(fig: Figure3D, traces: Any) -> Figure3D:
    if isinstance(traces, dict):
        for key, val in traces.items():
            if key == "payload" or not isinstance(val, dict):
                continue
            if val.get("trace_type") or val.get("type"):
                fig.add_trace(val)
        return fig
    if isinstance(traces, (list, tuple)):
        for t in traces:
            add_traces(fig, t)
    return fig

def plot_surface_dgrid(dgrid, isovalue, ase_atoms, region_labels=None, select="all", opacity=1.0,
                       cscale=None, color="blue", stepsize=1, stylize=False, random_seed=0,
                       hover_label="Region :", taubin_smooth=False, taubin_iterations=3,
                       taubin_lambda=0.2175, taubin_nu=-0.2175, **kwargs_mesh):
    del select, color, random_seed
    name = kwargs_mesh.pop("name", f"{isovalue} A isosurface")
    trace_id = kwargs_mesh.pop("trace_id", name.lower().replace(" ", "_"))
    regions = region_labels if region_labels is not None else np.zeros_like(dgrid, dtype=np.int32)
    style = _style_from_plotly_kwargs(opacity=opacity, stylize=stylize, name=name, **kwargs_mesh)
    return make_surface_trace(np.asarray(dgrid), np.asarray(regions), ase_atoms, isovalue=float(isovalue),
        step_size=int(stepsize), style=style, name=name, trace_id=trace_id, cscale=cscale,
        hover_prefix=hover_label.rstrip(" :"), taubin_smooth=bool(taubin_smooth),
        taubin_iterations=int(taubin_iterations), taubin_lambda=float(taubin_lambda), taubin_nu=float(taubin_nu))

def plot_cell(ase_atoms, color="purple", name="cell", **kwargs):
    return make_cell_trace(ase_atoms, style=_style_from_plotly_kwargs(color=color, name=name, **kwargs), name=name)

def draw_framework(ase_atoms, bonds="auto", representation="scatter", atom_scale=0.65,
                   bond_color=None, bond_linewidth=8.0, name="Framework", **kwargs):
    del kwargs
    if representation in ("scatter", "points"):
        return make_framework_scatter_trace(ase_atoms, bonds=bonds, name=name, atom_scale=atom_scale,
            bond_color=bond_color, bond_linewidth=bond_linewidth)
    from .movies import make_framework_trace
    return {"mesh": make_framework_trace(ase_atoms, bonds=bonds, name=name)}

def plot_maxima(seg_or_dgrid, atoms_or_maxima=None, ase_atoms=None, *, seg=None, hide_periodic=False,
                name="Maxima", cscale=None, **kwargs):
    del kwargs
    if seg is not None or isinstance(seg_or_dgrid, dict):
        seg_dict = prepare_seg_threejs(seg or seg_or_dgrid)
        atoms = atoms_or_maxima
    else:
        return plot_maxima_from_grid(seg_or_dgrid, atoms_or_maxima, ase_atoms, name=name, cscale=cscale, hide_periodic=hide_periodic)
    bundle = make_maxima_trace(seg_dict, atoms, name=name, cscale=cscale)
    if hide_periodic:
        bundle["outside"]["visible"] = False
        bundle["pbc"]["visible"] = False
    return bundle

def plot_maxima_from_grid(dgrid, maxima, ase_atoms, *, regions=None, name="Maxima", cscale=None, hide_periodic=False, **kwargs):
    del kwargs
    dgrid = np.asarray(dgrid)
    maxima = np.asarray(maxima, dtype=float)
    if regions is None:
        regions = np.zeros(dgrid.shape, dtype=np.int32)
        for i, w in enumerate(np.mod(maxima.astype(int), np.asarray(dgrid.shape, dtype=int))):
            regions[tuple(w)] = (i % 20) + 1
    seg_dict = {"dgrid": dgrid, "regions": np.asarray(regions), "all_maxima": maxima,
                "connections": [], "maxima_pbc": np.asarray([]), "maxima_info": np.asarray([]), "maxima_pbc_dict": {}}
    return plot_maxima(seg_dict, ase_atoms, name=name, cscale=cscale, hide_periodic=hide_periodic)

def plot_connections(seg, atoms, hide_periodic=True, **kwargs):
    del kwargs
    bundle = make_connections_trace(prepare_seg_threejs(seg), atoms)
    if hide_periodic:
        bundle["periodic_lines"]["visible"] = False
        bundle["periodic_nodes"]["visible"] = False
    return bundle

def plot_connections_as_it_is(connections, maxima, shape, ase_atoms, plot_periodic=True, plot_disconnected=False, **kwargs):
    del plot_disconnected
    seg = seg_from_connections(connections, maxima, shape)
    payload = _serialize_connection_traces(seg, ase_atoms)
    name = kwargs.pop("name", "Connections")
    trace_id = kwargs.pop("trace_id", "connections")
    style_lines = kwargs.pop("style", "internal_connections")
    def _line(suffix, pos_key, hover_key, periodic):
        positions = payload[pos_key]
        return _make_trace(trace_id=f"{trace_id}_{suffix}", trace_type="line_segments",
            name=f"{name} ({'periodic' if periodic else 'internal'})",
            data={"positions": positions, "n_segments": len(positions)//6},
            style="periodic_connections" if periodic else style_lines, hover=payload[hover_key],
            visible=plot_periodic if periodic else True)
    return {"internal_lines": _line("internal","internal_segments","internal_hover",False),
            "periodic_lines": _line("periodic","periodic_segments","periodic_hover",True),
            "internal_nodes": _make_trace(trace_id=f"{trace_id}_internal_nodes", trace_type="points",
                name=f"{name} internal nodes", data={"positions": payload["internal_points"],
                "n_points": len(payload["internal_points"])//3}, style="internal_nodes",
                hover=payload["internal_points_hover"], visible=True),
            "periodic_nodes": _make_trace(trace_id=f"{trace_id}_periodic_nodes", trace_type="points",
                name=f"{name} periodic nodes", data={"positions": payload["periodic_points"],
                "n_points": len(payload["periodic_points"])//3}, style="periodic_nodes",
                hover=payload["periodic_points_hover"], visible=plot_periodic), "payload": payload}

def plot_windows(connections_or_seg, name="Windows", marker_color=None, marker_size=10, marker_opacity=0.8, **kwargs):
    del marker_size
    if isinstance(connections_or_seg, dict):
        seg = prepare_seg_threejs(connections_or_seg)
    else:
        shape = kwargs.pop("grid_shape", None)
        maxima = kwargs.pop("maxima", [])
        if shape is None:
            raise ValueError("plot_windows(connections, ...) needs grid_shape= and maxima=")
        seg = seg_from_connections(connections_or_seg, maxima, shape)
    bundle = make_windows_trace(seg, name=name)
    if marker_color is not None:
        rgb = _hex_to_rgb01(marker_color)
        for key in ("internal", "periodic"):
            bundle[key].setdefault("style", {})["color"] = rgb
            bundle[key]["style"]["opacity"] = float(marker_opacity)
    return bundle


# ---------------------------------------------------------------------------
# Figure3JS — plotly-style chaining API for Three.js scenes
# ---------------------------------------------------------------------------

class Figure3JS:
    """Plotly-like granular builder for Three.js MOF visualisations.

    Mirror of the plotly trace API from ``mgr.visualize`` (``plot_surface_dgrid``,
    ``plot_cell``, ``plot_maxima``, ``plot_connections_as_it_is``, ``plot_windows``,
    ``draw_framework``) but targeting the Three.js backend (``MovieFigure``).

    All ``plot_*`` / ``draw_*`` methods return ``self`` so they can be chained::

        fig = (
            Figure3JS()
            .plot_framework(atoms)
            .plot_cell(atoms, color="orange")
            .plot_surface(dgrid, regions, atoms, name="Regions", opacity=1.0)
            .plot_surface(dgrid, regions_pbc, atoms, name="Regions PBC",
                          trace_id="surface_pbc", opacity=0.35)
            .plot_maxima(dgrid, all_maxima, atoms, name="Maxima", cscale=cscale)
            .plot_connections(connections, all_maxima, dgrid.shape, atoms)
            .plot_windows(connections, atoms=atoms, name="Windows",
                          marker_color="yellow")
            .to_html("output.html")
        )

    The underlying ``MovieFigure`` is always accessible via ``.figure``.

    Parameters
    ----------
    layout : "dark" | "white"
        Scene background preset.
    stylize : bool
        Apply PBR material tweaks (roughness / metalness) to mesh traces.
    """

    def __init__(self, layout: str = "dark", stylize: bool = True) -> None:
        self._fig = MovieFigure()
        self._stylize = stylize
        self._layout = layout
        if layout == "white":
            self._fig.set_appearance(background="#ffffff", scene_bg=0xFFFFFF)
        else:
            self._fig.set_appearance(
                background="radial-gradient(circle at center, #7a7a7a 0%, #111111 100%)",
                scene_bg=0x111111,
            )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _build_seg(
        dgrid: "np.ndarray | None" = None,
        regions: "np.ndarray | None" = None,
        all_maxima: "np.ndarray | None" = None,
        maxima_pbc: "np.ndarray | None" = None,
        maxima_info: "np.ndarray | None" = None,
        maxima_pbc_dict: "dict | None" = None,
        connections: "list | None" = None,
    ) -> "dict[str, Any]":
        """Build a minimal seg dict from raw user-facing arrays."""
        _dgrid = np.asarray(dgrid) if dgrid is not None else np.zeros((2, 2, 2), dtype=np.float32)
        _regions = np.asarray(regions) if regions is not None else np.zeros(_dgrid.shape, dtype=np.int32)
        return {
            "dgrid": _dgrid,
            "regions": _regions,
            "all_maxima": np.asarray(all_maxima, dtype=float) if all_maxima is not None else np.zeros((0, 3), dtype=float),
            "maxima_pbc": np.asarray(maxima_pbc, dtype=int) if maxima_pbc is not None else np.zeros((0, 3), dtype=int),
            "maxima_info": np.asarray(maxima_info, dtype=int) if maxima_info is not None else np.zeros((0, 7), dtype=int),
            "maxima_pbc_dict": maxima_pbc_dict if maxima_pbc_dict is not None else {},
            "connections": list(connections) if connections is not None else [],
        }

    def _maybe_stylize(self) -> None:
        if self._stylize:
            self._fig.update_traces(
                {"roughness": 0.05, "metalness": 0.3}, selector={"type": "mesh"}
            )

    def _apply_cscale(self, cscale: "list | None") -> "list[list[float]] | None":
        """Parse a plotly/raw colorscale and register it on the figure; return float form.

        Accepts:
        - list of plotly hex/rgb strings  → parsed to [[r,g,b], …] in [0,1]
        - list of [r,g,b] floats in [0,1] → used as-is
        - list of [r,g,b] ints in [0,255] → normalised to [0,1]
        """
        if cscale is None:
            return None
        first = cscale[0] if cscale else None
        if first is None:
            return None
        if isinstance(first, str):
            # plotly hex/rgb strings path — let MovieFigure.parse_plotly_colors handle it
            try:
                parsed = MovieFigure.parse_plotly_colors(cscale)
            except Exception:
                parsed = [[float(v)/255 for v in _hex_to_rgb01(c)] if isinstance(c, str) else [float(v) for v in c]
                          for c in cscale]
        else:
            # Already numeric: [[r,g,b], …] — normalise if values > 1 (assumed 0-255)
            arr = np.asarray(first, dtype=float)
            scale = 255.0 if arr.max() > 1.0 else 1.0
            parsed = [[float(v) / scale for v in rgb] for rgb in cscale]
        # Store directly in layout (bypasses set_colorscale which re-parses strings)
        if "colorscale" not in self._fig.layout:
            self._fig.layout["colorscale"] = {}
        self._fig.layout["colorscale"]["surface"] = parsed
        return parsed

    # ------------------------------------------------------------------
    # Framework
    # ------------------------------------------------------------------

    def plot_framework(
        self,
        atoms: "Atoms",
        bonds: "str | list | None" = "auto",
        representation: str = "ball_and_stick",
        name: "str | None" = None,
        **style: Any,
    ) -> "Figure3JS":
        """Add framework mesh (ball-and-stick, stick, or space-filling).

        Mirrors ``mgr.draw_framework()``.
        """
        kw: "dict[str, Any]" = {}
        if name is not None:
            kw["name"] = name
        self._fig.add_framework(atoms, bonds=bonds, representation=representation, **kw, **style)
        self._maybe_stylize()
        return self

    # Alias matching the plotly API name
    draw_framework = plot_framework  # type: ignore[assignment]

    # ------------------------------------------------------------------
    # Unit-cell wireframe
    # ------------------------------------------------------------------

    def plot_cell(
        self,
        atoms: "Atoms",
        color: "str | None" = None,
        name: str = "cell",
        **style: Any,
    ) -> "Figure3JS":
        """Add unit-cell wireframe.  Mirrors ``mgr.plot_cell()``."""
        if color is not None:
            style = {"color": _hex_to_rgb01(color), **style}
        self._fig.add_cell(atoms, name=name, **style)
        return self

    # ------------------------------------------------------------------
    # Isosurface / segmentation
    # ------------------------------------------------------------------

    def plot_surface(
        self,
        dgrid: "np.ndarray",
        region_labels: "np.ndarray",
        atoms: "Atoms",
        isovalue: float = 1.0,
        name: str = "Regions",
        opacity: float = 1.0,
        trace_id: "str | None" = None,
        cscale: "list | None" = None,
        step_size: int = 1,
        stylize: bool = False,
        taubin_smooth: bool = True,
        taubin_iterations: int = 3,
        taubin_lambda: float = 0.2175,
        taubin_nu: float = -0.2175,
        hover_label: str = "Region",
        **style: Any,
    ) -> "Figure3JS":
        """Add segmentation isosurface mesh.

        Mirrors ``mgr.plot_surface_dgrid()``.  Unlike the functional version
        this method *adds* the trace directly to the figure and returns ``self``.
        """
        if trace_id is None:
            trace_id = name.lower().replace(" ", "_").replace("(", "").replace(")", "").replace("-", "_")
        cs = self._apply_cscale(cscale)
        _style = _style_from_plotly_kwargs(opacity=opacity, stylize=stylize, name=name, **style)
        trace = make_surface_trace(
            np.asarray(dgrid),
            np.asarray(region_labels),
            atoms,
            isovalue=float(isovalue),
            step_size=int(step_size),
            style=_style,
            name=name,
            trace_id=trace_id,
            cscale=cs,
            hover_prefix=hover_label.rstrip(" :"),
            taubin_smooth=bool(taubin_smooth),
            taubin_iterations=int(taubin_iterations),
            taubin_lambda=float(taubin_lambda),
            taubin_nu=float(taubin_nu),
        )
        self._fig.add_trace(trace)
        self._fig.update_traces({"opacity": float(opacity)}, selector={"id": trace_id})
        # Non-default surfaces get a legend entry
        if trace_id != "surface":
            self._fig.groups[trace_id] = {"traces": [trace_id], "label": name}
        if self._stylize:
            self._maybe_stylize()
        return self

    # Alias matching the functional API
    def plot_surface_dgrid(
        self,
        dgrid: "np.ndarray",
        isovalue: float,
        ase_atoms: "Atoms",
        region_labels: "np.ndarray | None" = None,
        opacity: float = 1.0,
        cscale: "list | None" = None,
        stepsize: int = 1,
        stylize: bool = False,
        taubin_smooth: bool = True,
        taubin_iterations: int = 3,
        taubin_lambda: float = 0.2175,
        taubin_nu: float = -0.2175,
        **kwargs: Any,
    ) -> "Figure3JS":
        """Functional-API alias for ``plot_surface``; matches ``mgr.plot_surface_dgrid()``."""
        name = kwargs.pop("name", f"{isovalue} A isosurface")
        trace_id = kwargs.pop("trace_id", None)
        regions = region_labels if region_labels is not None else np.zeros_like(np.asarray(dgrid), dtype=np.int32)
        return self.plot_surface(
            dgrid, regions, ase_atoms,
            isovalue=float(isovalue),
            name=name,
            opacity=float(opacity),
            trace_id=trace_id,
            cscale=cscale,
            step_size=int(stepsize),
            stylize=stylize,
            taubin_smooth=taubin_smooth,
            taubin_iterations=taubin_iterations,
            taubin_lambda=taubin_lambda,
            taubin_nu=taubin_nu,
            **kwargs,
        )

    # ------------------------------------------------------------------
    # Maxima (pore centres)
    # ------------------------------------------------------------------

    def plot_maxima(
        self,
        dgrid_or_seg: "np.ndarray | dict[str, Any]",
        maxima_or_atoms: "np.ndarray | Atoms | None" = None,
        ase_atoms: "Atoms | None" = None,
        *,
        regions: "np.ndarray | None" = None,
        maxima_pbc: "np.ndarray | None" = None,
        maxima_info: "np.ndarray | None" = None,
        maxima_pbc_dict: "dict | None" = None,
        connections: "list | None" = None,
        name: str = "Maxima",
        trace_id: str = "maxima",
        cscale: "list | None" = None,
        hide_periodic: bool = True,
        **style: Any,
    ) -> "Figure3JS":
        """Add pore-centre maxima scatter points.

        Accepts either:
        - ``(seg_dict, atoms)``  — seg-dict path, matches ``mgr.plot_maxima(seg, atoms)``
        - ``(dgrid, maxima_arr, atoms)`` — raw arrays, matches ``mgr.plot_maxima(dgrid, maxima, atoms)``

        Mirrors ``mgr.plot_maxima()`` / ``mgr.plot_maxima_from_grid()``.
        """
        cs = self._apply_cscale(cscale)
        if isinstance(dgrid_or_seg, dict):
            seg = prepare_seg_threejs(dgrid_or_seg)
            atoms = maxima_or_atoms
        else:
            # Raw arrays path
            dgrid = np.asarray(dgrid_or_seg)
            maxima = np.asarray(maxima_or_atoms, dtype=float)
            atoms = ase_atoms
            if regions is None:
                _regions = np.zeros(dgrid.shape, dtype=np.int32)
                for _i, _w in enumerate(np.mod(maxima.astype(int), np.asarray(dgrid.shape, dtype=int))):
                    _regions[tuple(_w)] = (_i % 20) + 1
            else:
                _regions = np.asarray(regions)
            seg = self._build_seg(
                dgrid=dgrid, regions=_regions, all_maxima=maxima,
                maxima_pbc=maxima_pbc, maxima_info=maxima_info,
                maxima_pbc_dict=maxima_pbc_dict, connections=connections,
            )
        bundle = make_maxima_trace(seg, atoms, name=name, trace_id=trace_id, cscale=cs)
        if hide_periodic:
            bundle["outside"]["visible"] = False
            bundle["pbc"]["visible"] = False
        self._fig.add_trace(bundle["internal"])
        self._fig.add_trace(bundle["outside"])
        self._fig.add_trace(bundle["pbc"])
        # Legend groups
        self._fig.groups[f"{trace_id}_internal"] = {
            "traces": [f"{trace_id}_internal"], "label": f"{name} (internal)",
        }
        if not hide_periodic:
            self._fig.groups[f"{trace_id}_pbc"] = {
                "traces": [f"{trace_id}_pbc"], "label": f"{name} (PBC reps)",
            }
        return self

    # Functional-API alias
    def plot_maxima_from_grid(
        self,
        dgrid: "np.ndarray",
        maxima: "np.ndarray",
        ase_atoms: "Atoms",
        *,
        regions: "np.ndarray | None" = None,
        name: str = "Maxima",
        trace_id: str = "maxima",
        cscale: "list | None" = None,
        hide_periodic: bool = False,
        **kwargs: Any,
    ) -> "Figure3JS":
        """Raw-arrays alias; mirrors ``mgr.plot_maxima_from_grid()``."""
        return self.plot_maxima(
            dgrid, maxima, ase_atoms,
            regions=regions, name=name, trace_id=trace_id,
            cscale=cscale, hide_periodic=hide_periodic, **kwargs,
        )

    # ------------------------------------------------------------------
    # Connections (pore graph)
    # ------------------------------------------------------------------

    def plot_connections(
        self,
        connections: list,
        all_maxima: "np.ndarray",
        shape: "tuple | np.ndarray",
        atoms: "Atoms",
        *,
        dgrid: "np.ndarray | None" = None,
        regions: "np.ndarray | None" = None,
        hide_periodic: bool = True,
        plot_periodic: bool = False,
        name: str = "Connections",
        trace_id: str = "connections",
        **style: Any,
    ) -> "Figure3JS":
        """Add pore-graph connection lines and node markers.

        Mirrors ``mgr.plot_connections_as_it_is()``.
        """
        seg = seg_from_connections(connections, all_maxima, shape, dgrid=dgrid)
        bundle = plot_connections_as_it_is(
            connections, all_maxima, shape, atoms,
            plot_periodic=plot_periodic, name=name, trace_id=trace_id, **style,
        )
        add_traces(self._fig, bundle)
        return self

    # Alias for full sig-match with the plotly API
    def plot_connections_as_it_is(
        self,
        connections: list,
        all_maxima: "np.ndarray",
        shape: "tuple | np.ndarray",
        ase_atoms: "Atoms",
        *,
        plot_periodic: bool = True,
        plot_disconnected: bool = False,
        name: str = "Connections",
        trace_id: str = "connections",
        cscale: "list | None" = None,
        **style: Any,
    ) -> "Figure3JS":
        """Exact-signature alias matching ``mgr.plot_connections_as_it_is()``."""
        bundle = plot_connections_as_it_is(
            connections, all_maxima, shape, ase_atoms,
            plot_periodic=plot_periodic, plot_disconnected=plot_disconnected,
            name=name, trace_id=trace_id, **style,
        )
        add_traces(self._fig, bundle)
        return self

    # ------------------------------------------------------------------
    # Windows
    # ------------------------------------------------------------------

    def plot_windows(
        self,
        connections: "list | dict",
        atoms: "Atoms | None" = None,
        name: str = "Windows",
        trace_id: str = "windows",
        marker_color: "str | None" = None,
        marker_size: float = 10,
        marker_opacity: float = 0.8,
        hide_periodic: bool = True,
        spacing: "float | None" = None,
        dgrid: "np.ndarray | None" = None,
        maxima: "np.ndarray | None" = None,
        grid_shape: "tuple | None" = None,
        **kwargs: Any,
    ) -> "Figure3JS":
        """Add window markers between connected pores.

        Mirrors ``mgr.plot_windows()``.  Accepts either a seg-dict or a raw
        connections list (requires ``grid_shape=`` and ``maxima=``).
        """
        bundle = plot_windows(
            connections,
            name=name,
            marker_color=marker_color,
            marker_size=marker_size,
            marker_opacity=marker_opacity,
            maxima=maxima if maxima is not None else [],
            grid_shape=grid_shape,
        )
        if hide_periodic:
            bundle.get("periodic", {})["visible"] = False
        add_traces(self._fig, bundle)
        return self

    # ------------------------------------------------------------------
    # Guest / adsorbate atoms
    # ------------------------------------------------------------------

    def plot_guests(
        self,
        atoms: "Atoms",
        color: "list | None" = None,
        **style: Any,
    ) -> "Figure3JS":
        """Add guest/adsorbate atoms as instanced spheres."""
        self._fig.add_guests(atoms, color=color, **style)
        return self

    # ------------------------------------------------------------------
    # Layout / camera / appearance
    # ------------------------------------------------------------------

    def set_camera(
        self,
        eye: "tuple | None" = None,
        center: "tuple | None" = None,
        up: "tuple | None" = None,
        projection: str = "perspective",
    ) -> "Figure3JS":
        """Set camera position and projection type."""
        self._fig.set_camera(eye=eye, center=center, up=up, projection=projection)
        return self

    def set_appearance(self, **kwargs: Any) -> "Figure3JS":
        """Set scene appearance (background colour, theme, …)."""
        self._fig.set_appearance(**kwargs)
        return self

    def update_traces(
        self,
        patch: "dict[str, Any]",
        selector: "dict[str, Any] | None" = None,
    ) -> "Figure3JS":
        """Apply a style patch to matching traces (opacity, color, visible, …)."""
        self._fig.update_traces(patch, selector=selector)
        return self

    def add_group(
        self,
        group_id: str,
        traces: "list[str]",
        label: str,
    ) -> "Figure3JS":
        """Manually register a legend group."""
        self._fig.groups[group_id] = {"traces": traces, "label": label}
        return self

    def add_trace(self, trace: "dict[str, Any]") -> "Figure3JS":
        """Add a raw trace dict directly (escape hatch for advanced use)."""
        self._fig.add_trace(trace)
        return self

    def add_traces(self, traces: Any) -> "Figure3JS":
        """Add a trace dict or bundle (dict-of-traces / list of traces)."""
        add_traces(self._fig, traces)
        return self

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def to_html(
        self,
        path: "str | 'Path' | None" = None,
        **kwargs: Any,
    ) -> "str | 'Path'":
        """Export to standalone HTML.

        Parameters
        ----------
        path : str or Path, optional
            Output file path.  If *None* the HTML string is returned instead.
        **kwargs
            Forwarded to ``MovieFigure.to_html()``.

        Returns
        -------
        str
            HTML string when *path* is None.
        Path
            Resolved output path when *path* is given.
        """
        if path is None:
            import tempfile, os
            tmp = tempfile.NamedTemporaryFile(suffix=".html", delete=False)
            tmp.close()
            self._fig.to_html(tmp.name, **kwargs)
            with open(tmp.name, "r", encoding="utf-8") as f:
                html = f.read()
            os.unlink(tmp.name)
            return html
        return self._fig.to_html(path, **kwargs)

    def show(self, renderer: str = "browser") -> None:
        """Display figure in a browser tab or Jupyter notebook."""
        self._fig.show(renderer=renderer)

    def _repr_html_(self) -> str:  # Jupyter integration
        return self._fig._repr_html_()

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def figure(self) -> "MovieFigure":
        """The underlying ``MovieFigure`` for advanced / low-level access."""
        return self._fig

    @property
    def data(self) -> list:
        return self._fig.data

    @property
    def groups(self) -> dict:
        return self._fig.groups
