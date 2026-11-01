"""Widest-path bottleneck re-partition for leaky PBC marker pairs.

Phase 1 prototype: fix rare wrap-around pocket segmentation where standard
watershed places the periodic self-boundary at the cell face instead of the
internal constriction.
"""

from __future__ import annotations

import heapq
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

Offset = Tuple[int, int, int]


def neighbor_offsets_26() -> List[Offset]:
    """Return the 26 face/edge/corner neighbor offsets for a 3D voxel grid."""
    return [
        (dx, dy, dz)
        for dx in (-1, 0, 1)
        for dy in (-1, 0, 1)
        for dz in (-1, 0, 1)
        if not (dx == 0 and dy == 0 and dz == 0)
    ]


def marker_padded_position(
    row: Sequence[int],
    dgrid_shape: Tuple[int, int, int],
    padding: int,
) -> Tuple[int, int, int]:
    """Padded-grid index for one ``maxima_info`` row."""
    orig = np.array(dgrid_shape, dtype=int)
    pos = np.asarray(row[1:4], dtype=int) + int(padding) + np.asarray(row[4:7], dtype=int) * orig
    return int(pos[0]), int(pos[1]), int(pos[2])


def find_all_pbc_wrap_marker_pairs(
    maxima_info: np.ndarray,
    pbc_groups: Optional[Sequence] = None,
) -> List[Tuple[int, int, Offset]]:
    """All in-cell + single-axis wrapped label pairs at each maxima position.

    Pairs are found by shared unit-cell image coordinates ``row[1:4]``, not only
    ``pbc_groups`` membership (a-axis wraps like label 17 may sit outside the group).
    """
    info = np.asarray(maxima_info, dtype=int)
    by_pos: Dict[Tuple[int, int, int], List[np.ndarray]] = {}
    for row in info:
        key = (int(row[1]), int(row[2]), int(row[3]))
        by_pos.setdefault(key, []).append(row)

    pairs: List[Tuple[int, int, Offset]] = []
    for rows in by_pos.values():
        in_cell = [r for r in rows if tuple(r[4:7]) == (0, 0, 0)]
        if not in_cell:
            continue
        r0 = in_cell[0]
        for rw in rows:
            off = tuple(int(x) for x in rw[4:7])
            if off == (0, 0, 0):
                continue
            if sum(1 for x in off if x != 0) != 1:
                continue
            pairs.append((int(r0[0]), int(rw[0]), off))
    return pairs


def find_pbc_wrap_marker_pair(
    maxima_info: np.ndarray,
    pbc_groups: Sequence,
    labels: Optional[np.ndarray] = None,
    dgrid: Optional[np.ndarray] = None,
    mask: Optional[np.ndarray] = None,
    boundary_score_fn=None,
) -> Optional[Tuple[int, int, Tuple[int, int, int], Tuple[int, int, int], Offset]]:
    """Find the leakiest in-cell + single-axis wrapped marker pair.

    When ``labels``, ``dgrid``, ``mask``, and ``boundary_score_fn`` are supplied,
    picks the pair with the largest boundary-score (Phase 2 leak detector).
    Otherwise returns the first single-axis pair found.

    Returns ``(label_a, label_b, pos_a, pos_b, wrap_offset)`` or ``None``.
    """
    info = np.asarray(maxima_info, dtype=int)
    _ = pbc_groups  # retained for API compatibility; pairs use image position
    candidates = find_all_pbc_wrap_marker_pairs(info)
    if not candidates:
        return None

    if (
        labels is not None
        and dgrid is not None
        and mask is not None
        and boundary_score_fn is not None
    ):
        best = None
        best_score = -1.0
        for la, lb, off in candidates:
            score = float(boundary_score_fn(labels, dgrid, mask, la, lb))
            if score > best_score + 1e-6:
                best_score = score
                best = (la, lb, off)
            elif best is not None and abs(score - best_score) <= 1e-6:
                # Tie-break: prefer +axis primary wrap (e.g. (1,0,0) over (-1,0,0))
                if off > best[2]:
                    best = (la, lb, off)
        if best is None:
            return None
        la, lb, off = best
    else:
        la, lb, off = candidates[0]

    row_a = info[info[:, 0] == la][0]
    row_b = info[info[:, 0] == lb][0]
    return (
        la,
        lb,
        (int(row_a[1]), int(row_a[2]), int(row_a[3])),
        (int(row_b[1]), int(row_b[2]), int(row_b[3])),
        off,
    )


def widest_path_bottleneck_from_seed(
    dgrid: np.ndarray,
    mask: np.ndarray,
    seed: Tuple[int, int, int],
    offsets: Optional[Sequence[Offset]] = None,
) -> np.ndarray:
    """Max-bottleneck distances from ``seed`` on the void graph (26-NN).

    ``bottleneck[v] = max over paths seed→v of min(dgrid along path)``.
    Unreachable voxels remain ``-1``.
    """
    if offsets is None:
        offsets = neighbor_offsets_26()

    shape = dgrid.shape
    nx, ny, nz = shape
    si, sj, sk = seed
    if not (0 <= si < nx and 0 <= sj < ny and 0 <= sk < nz):
        raise ValueError(f"seed {seed} out of bounds for shape {shape}")
    if not mask[si, sj, sk]:
        raise ValueError(f"seed {seed} not in void mask")

    best = np.full(shape, -1.0, dtype=np.float32)
    start = float(dgrid[si, sj, sk])
    best[si, sj, sk] = start
    heap: List[Tuple[float, int, int, int]] = [(-start, si, sj, sk)]

    while heap:
        neg_b, i, j, k = heapq.heappop(heap)
        b = -neg_b
        if b < float(best[i, j, k]) - 1e-7:
            continue

        for dx, dy, dz in offsets:
            ni, nj, nk = i + dx, j + dy, k + dz
            if ni < 0 or nj < 0 or nk < 0 or ni >= nx or nj >= ny or nk >= nz:
                continue
            if not mask[ni, nj, nk]:
                continue
            cand = min(b, float(dgrid[ni, nj, nk]))
            if cand > float(best[ni, nj, nk]) + 1e-7:
                best[ni, nj, nk] = cand
                heapq.heappush(heap, (-cand, ni, nj, nk))

    return best


def repartition_pair_widest_path(
    labels: np.ndarray,
    dgrid: np.ndarray,
    mask: np.ndarray,
    label_a: int,
    label_b: int,
    seed_a: Tuple[int, int, int],
    seed_b: Tuple[int, int, int],
    offsets: Optional[Sequence[Offset]] = None,
) -> Tuple[np.ndarray, Dict[str, float]]:
    """Reassign disputed voxels between two seeds via widest-path partition."""
    if offsets is None:
        offsets = neighbor_offsets_26()

    out = np.asarray(labels, dtype=np.uint32, order="C").copy()
    b_a = widest_path_bottleneck_from_seed(dgrid, mask, seed_a, offsets)
    b_b = widest_path_bottleneck_from_seed(dgrid, mask, seed_b, offsets)

    reachable_a = b_a >= 0
    reachable_b = b_b >= 0
    either_label = (out == np.uint32(label_a)) | (out == np.uint32(label_b))
    disputed = mask & reachable_a & reachable_b & either_label

    n_disputed = int(np.count_nonzero(disputed))
    if n_disputed == 0:
        return out, {
            "disputed_voxels": 0.0,
            "reassigned_to_a": 0.0,
            "reassigned_to_b": 0.0,
            "unchanged": 0.0,
        }

    assign_a = disputed & (b_a > b_b)
    assign_b = disputed & (b_b > b_a)
    tie = disputed & (~assign_a) & (~assign_b)
    if np.any(tie):
        if label_a <= label_b:
            assign_a = assign_a | tie
        else:
            assign_b = assign_b | tie

    before_a = int(np.count_nonzero((out == label_a) & disputed))
    before_b = int(np.count_nonzero((out == label_b) & disputed))

    out[assign_a] = np.uint32(label_a)
    out[assign_b] = np.uint32(label_b)

    after_a = int(np.count_nonzero((out == label_a) & disputed))
    after_b = int(np.count_nonzero((out == label_b) & disputed))

    return out, {
        "disputed_voxels": float(n_disputed),
        "reassigned_to_a": float(after_a - before_a),
        "reassigned_to_b": float(after_b - before_b),
        "before_a_in_disputed": float(before_a),
        "before_b_in_disputed": float(before_b),
        "after_a_in_disputed": float(after_a),
        "after_b_in_disputed": float(after_b),
    }


def widest_path_bottleneck_between_seeds(
    dgrid: np.ndarray,
    mask: np.ndarray,
    seed_a: Tuple[int, int, int],
    seed_b: Tuple[int, int, int],
    offsets: Optional[Sequence[Offset]] = None,
) -> float:
    """Bottleneck of the widest path between two seeds (symmetric max-min)."""
    b_a = widest_path_bottleneck_from_seed(dgrid, mask, seed_a, offsets)
    si, sj, sk = seed_b
    return float(b_a[si, sj, sk])


def apply_widest_path_fix_to_pair(
    labels: np.ndarray,
    dgrid: np.ndarray,
    mask: np.ndarray,
    maxima_info: np.ndarray,
    pbc_groups: Sequence,
    dgrid_shape: Tuple[int, int, int],
    padding: int,
    boundary_score_fn=None,
) -> Tuple[np.ndarray, Optional[Dict]]:
    """Locate PBC wrap pair and re-partition; return labels unchanged if none found."""
    pair = find_pbc_wrap_marker_pair(
        maxima_info,
        pbc_groups,
        labels=labels,
        dgrid=dgrid,
        mask=mask,
        boundary_score_fn=boundary_score_fn,
    )
    if pair is None:
        return labels, None

    label_a, label_b, _img_a, _img_b, wrap_off = pair
    seed_a = marker_padded_position(
        maxima_info[maxima_info[:, 0] == label_a][0],
        dgrid_shape,
        padding,
    )
    seed_b = marker_padded_position(
        maxima_info[maxima_info[:, 0] == label_b][0],
        dgrid_shape,
        padding,
    )

    fixed, stats = repartition_pair_widest_path(
        labels=labels,
        dgrid=dgrid,
        mask=mask,
        label_a=label_a,
        label_b=label_b,
        seed_a=seed_a,
        seed_b=seed_b,
    )
    wr = widest_path_bottleneck_between_seeds(dgrid, mask, seed_a, seed_b)
    meta = {
        "label_a": label_a,
        "label_b": label_b,
        "seed_a": seed_a,
        "seed_b": seed_b,
        "wrap_offset": wrap_off,
        "widest_path_wr": wr,
        **stats,
    }
    return fixed, meta
