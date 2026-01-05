"""
Experimental window refinement utilities.

These helpers operate on the `results` dict returned by
`mofography.do_periodic_segmentation_rust` and produce refined window
locations/radii using different geometric heuristics.

The functions do **not** mutate the original results; instead they return
new connection arrays that mirror the structure of the original
`results["connections"]` entries:

    [array([reg_i, reg_j]), connected_flag, "internal"/"periodic",
     numpy.ndarray(window_centers), numpy.ndarray(window_radii)]
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import scipy.ndimage as ndi


@dataclass
class WindowRefinementOptions:
    """Common options for refinement routines."""

    search_radius_vox: int = 6
    require_unique_regions: bool = True
    radius_steps: int = 6


def refine_windows_via_morph_opening(
    results: Dict,
    *,
    cell_matrix: Optional[np.ndarray] = None,
    options: WindowRefinementOptions | None = None,
) -> List[np.ndarray]:
    """
    Refine window centers/radii by searching for the smallest distance-field
    threshold that disconnects the two regions (local morphological opening).
    """

    opts = options or WindowRefinementOptions()
    ctx = _RefinementContext.from_results(
        results, opts.require_unique_regions, cell_matrix=cell_matrix
    )
    refined = []
    for conn in ctx.connections:
        regs = conn[0].astype(int)
        if not ctx.should_refine(regs):
            refined.append(conn.copy())
            continue
        centers = np.asarray(conn[3])
        radii = np.asarray(conn[4])
        new_centers: List[np.ndarray] = []
        new_radii: List[float] = []
        for center_cart, radius in zip(centers, radii):
            refined_window = _refine_window_morph(
                ctx,
                center_cart,
                tuple(regs),
                search_radius=opts.search_radius_vox,
                n_steps=opts.radius_steps,
            )
            if refined_window is None:
                new_centers.append(center_cart)
                new_radii.append(radius)
            else:
                new_centers.append(refined_window.center)
                new_radii.append(refined_window.radius)
        refined.append(
            np.array(
                [
                    conn[0],
                    conn[1],
                    conn[2],
                    np.asarray(new_centers, dtype=float),
                    np.asarray(new_radii, dtype=float),
                ],
                dtype=object,
            )
        )
    return refined


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


@dataclass
class _RefinedWindow:
    center: np.ndarray
    radius: float


@dataclass
class _RefinementContext:
    dgrid: np.ndarray
    dgrid_padded: np.ndarray
    regions_padded: np.ndarray
    padding: int
    cell: np.ndarray
    _cell_inv: np.ndarray
    connections: Sequence[np.ndarray]
    unique_regions: Optional[np.ndarray]

    @classmethod
    def from_results(
        cls,
        results: Dict,
        require_unique_regions: bool,
        *,
        cell_matrix: Optional[np.ndarray] = None,
    ) -> "_RefinementContext":
        if "dgrid" not in results:
            raise ValueError(
                "results dictionary does not contain 'dgrid'. "
                "Re-run periodic segmentation with return_dgrid=True."
            )
        if "regions_padded" not in results:
            raise ValueError(
                "results dictionary does not contain 'regions_padded'. "
                "Re-run periodic segmentation with return_padded_regions=True."
            )
        dgrid = np.asarray(results["dgrid"])
        regions_padded = np.asarray(results["regions_padded"])
        connections = results.get("connections")
        if connections is None:
            raise ValueError("results dictionary does not contain 'connections'.")
        padding = int(results["input_params"]["padding_thickness"])
        if cell_matrix is not None:
            cell = np.asarray(cell_matrix, dtype=float)
        else:
            cell = np.asarray(results["input_params"].get("cell_matrix", []), dtype=float)
            if cell.size == 0:
                raise ValueError(
                    "cell_matrix not provided. Pass `cell_matrix=mat_atoms.cell.array` "
                    "when calling the refinement helper."
                )
        unique_regions = None
        if require_unique_regions:
            unique_regions = np.asarray(results.get("unique_regions", []), dtype=int)
        return cls(
            dgrid=dgrid,
            dgrid_padded=np.pad(dgrid, padding, mode="wrap"),
            regions_padded=regions_padded,
            padding=padding,
            cell=cell,
            _cell_inv=np.linalg.inv(cell),
            connections=connections,
            unique_regions=unique_regions,
        )

    def should_refine(self, regs: Sequence[int]) -> bool:
        if self.unique_regions is None:
            return True
        return bool(np.any(np.isin(regs, self.unique_regions)))

    @property
    def inner_shape(self) -> Tuple[int, int, int]:
        return (
            self.dgrid.shape[0] - 2 * self.padding,
            self.dgrid.shape[1] - 2 * self.padding,
            self.dgrid.shape[2] - 2 * self.padding,
        )

    @property
    def cell_inv(self) -> np.ndarray:
        return self._cell_inv


def _refine_window_morph(
    ctx: _RefinementContext,
    center_cart: np.ndarray,
    regs: Tuple[int, int],
    *,
    search_radius: int,
    n_steps: int,
) -> Optional[_RefinedWindow]:
    slices, center_idx = _bounding_box(center_cart, ctx, search_radius)
    label_chunk = ctx.regions_padded[slices]
    mask_chunk = np.isin(label_chunk, regs)
    if not mask_chunk.any():
        return None

    dgrid_chunk = ctx.dgrid_padded[slices]
    values = dgrid_chunk[mask_chunk]
    if values.size == 0:
        return None

    low = 0.0
    high = float(values.max())
    threshold = 0.0
    for _ in range(max(1, n_steps)):
        mid = 0.5 * (low + high)
        if _regions_connected(mask_chunk, label_chunk, dgrid_chunk, regs, mid):
            low = mid
        else:
            high = mid
        threshold = mid

    candidate = mask_chunk & (dgrid_chunk >= threshold)
    reg1_mask = candidate & (label_chunk == regs[0])
    reg2_mask = candidate & (label_chunk == regs[1])
    if not reg1_mask.any() or not reg2_mask.any():
        return None

    reg1_center = _voxel_centroid(reg1_mask, slices)
    reg2_center = _voxel_centroid(reg2_mask, slices)
    center = 0.5 * (reg1_center + reg2_center)
    return _RefinedWindow(center=center, radius=threshold)


def _bounding_box(
    center_cart: np.ndarray,
    ctx: _RefinementContext,
    search_radius: int,
) -> Tuple[Tuple[slice, slice, slice], np.ndarray]:
    frac = ctx.cell_inv @ center_cart
    frac = frac % 1.0
    idx_float = np.array(frac) * np.array(ctx.inner_shape) + ctx.padding
    idx = np.round(idx_float).astype(int)
    slices = []
    center_idx = []
    for axis, dim in enumerate(ctx.dgrid_padded.shape):
        lo = max(idx[axis] - search_radius, 0)
        hi = min(idx[axis] + search_radius + 1, dim)
        slices.append(slice(lo, hi))
        center_idx.append(idx[axis] - lo)
    return (tuple(slices), np.array(center_idx, dtype=int))


def _regions_connected(
    mask_chunk: np.ndarray,
    label_chunk: np.ndarray,
    dgrid_chunk: np.ndarray,
    regs: Tuple[int, int],
    threshold: float,
) -> bool:
    candidate = mask_chunk & (dgrid_chunk >= threshold)
    seeds = np.argwhere(candidate & (label_chunk == regs[0]))
    if seeds.size == 0:
        return False
    visited = np.zeros(candidate.shape, dtype=bool)
    queue = [tuple(seeds[0])]
    visited[queue[0]] = True
    while queue:
        x, y, z = queue.pop()
        if label_chunk[(x, y, z)] == regs[1]:
            return True
        for nx, ny, nz in _six_neighbors((x, y, z)):
            if (
                0 <= nx < candidate.shape[0]
                and 0 <= ny < candidate.shape[1]
                and 0 <= nz < candidate.shape[2]
                and candidate[(nx, ny, nz)]
                and not visited[(nx, ny, nz)]
            ):
                visited[(nx, ny, nz)] = True
                queue.append((nx, ny, nz))
    return False


def _voxel_centroid(mask: np.ndarray, slices: Tuple[slice, slice, slice]) -> np.ndarray:
    coords = np.argwhere(mask)
    if coords.size == 0:
        return np.zeros(3)
    offsets = np.array([sl.start for sl in slices])
    global_coords = coords + offsets
    return global_coords.mean(axis=0)


def _voxel_to_cartesian(
    voxel_idx: np.ndarray,
    slices: Tuple[slice, slice, slice],
    ctx: _RefinementContext,
) -> np.ndarray:
    global_idx = voxel_idx + np.array([sl.start for sl in slices])
    frac = (global_idx - ctx.padding) / np.array(ctx.inner_shape)
    frac = np.clip(frac, 0.0, 1.0)
    return ctx.cell @ frac


def _six_neighbors(idx: Iterable[int]) -> List[np.ndarray]:
    x, y, z = idx
    return [
        np.array([x + 1, y, z]),
        np.array([x - 1, y, z]),
        np.array([x, y + 1, z]),
        np.array([x, y - 1, z]),
        np.array([x, y, z + 1]),
        np.array([x, y, z - 1]),
    ]


def _within_chunk(idx: np.ndarray, shape: Tuple[int, int, int]) -> bool:
    return (
        0 <= idx[0] < shape[0]
        and 0 <= idx[1] < shape[1]
        and 0 <= idx[2] < shape[2]
    )


