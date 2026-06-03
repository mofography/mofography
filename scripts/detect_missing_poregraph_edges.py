#!/usr/bin/env python
"""
Missing Edge Detector for Periodic Pore Graphs
================================================
Strategy: distance-shell completion.

For each pair of base-cell nodes that already have ≥1 edge, compute
Cartesian distances for all offsets in {-1,0,1}^3.  Cluster into
distance shells.  If an edge exists in a shell but others in the same
shell are absent, those are candidate missing edges.

Outputs results as .npz so segmentation only runs once.
"""

import sys
from pathlib import Path

import numpy as np

MOFOGRAPHY_ROOT = Path("/Users/arungopalan/repos/mofography")
sys.path.insert(0, str(MOFOGRAPHY_ROOT))

CIF = str(MOFOGRAPHY_ROOT / "examples" / "TestCIFs" / "NU-1000.cif")
CACHE = Path("/tmp/nu1000_missing_edges_cache.npz")

# ── helpers ───────────────────────────────────────────────────────────

def canonical_key(i, j, offset):
    """Return (min, max, canonical_offset) with i<j convention."""
    offset = tuple(int(x) for x in offset)
    if i < j:
        return (i, j, offset)
    elif i > j:
        return (j, i, tuple(-x for x in offset))
    else:
        # self-edge: canonicalize offset (positive lexicographic)
        if offset == (0, 0, 0):
            return None
        neg = tuple(-x for x in offset)
        if offset < neg:
            return (i, i, offset)
        else:
            return (i, i, neg)


def cart_distance(frac_i, frac_j, offset, cell):
    """Cartesian distance from node i at origin to node j at +offset."""
    d_frac = frac_j + np.array(offset, dtype=float) - frac_i
    d_cart = d_frac @ cell.T
    return np.linalg.norm(d_cart), d_cart


def cluster_shell(distances, tol=0.15):
    """
    Group (key, distance) pairs into shells by relative distance similarity.
    Returns list of shells, each shell is a list of (key, dist).
    Simple greedy clustering sorted by distance.
    """
    if not distances:
        return []
    items = sorted(distances, key=lambda x: x[1])
    shells = []
    current_shell = [items[0]]
    ref_dist = items[0][1]
    for key, d in items[1:]:
        if ref_dist > 0 and abs(d - ref_dist) / ref_dist < tol:
            current_shell.append((key, d))
        else:
            shells.append(current_shell)
            current_shell = [(key, d)]
            ref_dist = d
    shells.append(current_shell)
    return shells


# ── main ──────────────────────────────────────────────────────────────

def main():
    from ase.geometry import cellpar_to_cell
    from mofography.cif_io import read_cif_to_ase
    import mofography as mgr

    # ── Step 0: load or compute segmentation ──────────────────────────
    if CACHE.exists():
        print(f"Loading cached segmentation from {CACHE}")
        data = dict(np.load(CACHE, allow_pickle=True))
        maxima_info = data["maxima_info"]
        connections = list(data["connections"])
        dgrid_shape = tuple(data["dgrid_shape"])
        cellpar = data["cellpar"]
        existing_edges_arr = data["existing_edges"]
        frac_pos_arr = data["frac_pos"]
        pocket_radii_arr = data["pocket_radii"]
    else:
        print("Running segmentation …")
        atoms = read_cif_to_ase(CIF)
        atoms = mgr.update_with_covalent_bonds(atoms)
        mgr.compute_aromatic_cycles(atoms)

        dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=0.5)
        seg = mgr.do_periodic_segmentation_rust(
            dgrid, atoms, spacing=0.5, mask_thickness=0.05, h=0.2,
            threshold_abs=1.5, min_distance=1.5, extra_padding=5,
            compute_connections=True, return_dgrid=False)

        from tipmof.tipmof import collect_periodic_edges, make_directed_graph_from_segmentation_results
        G = make_directed_graph_from_segmentation_results(seg, dgrid=dgrid)
        edge_map, _ = collect_periodic_edges(seg["connections"], seg["maxima_info"])

        maxima_info = seg["maxima_info"]
        connections = seg["connections"]
        dgrid_shape = dgrid.shape
        cellpar = atoms.cell.cellpar()

        # Collect per-node data
        is_base = np.all(maxima_info[:, 4:7] == 0, axis=1)
        base_labels = sorted(int(maxima_info[i, 0]) for i in np.where(is_base)[0])

        frac_pos = {}
        pocket_radii = {}
        for node in G.nodes():
            frac_pos[node] = np.array(G.nodes[node]["frac_pos"], dtype=float)
            pocket_radii[node] = G.nodes[node].get("pocket_radius", 0.0)

        existing_edges = list(edge_map.keys())  # (i, j, offset) tuples

        # Save cache
        np.savez(
            CACHE,
            maxima_info=maxima_info,
            connections=np.array(connections, dtype=object),
            dgrid_shape=np.array(dgrid_shape),
            cellpar=cellpar,
            existing_edges=np.array(existing_edges, dtype=object),
            frac_pos=np.array([frac_pos[n] for n in sorted(frac_pos)], dtype=float),
            pocket_radii=np.array([pocket_radii[n] for n in sorted(pocket_radii)], dtype=float),
            base_labels=np.array(base_labels),
        )
        print(f"Cached segmentation to {CACHE}")

    # ── Step 1: reconstruct data structures ───────────────────────────
    cell = cellpar_to_cell(cellpar)
    existing_set = set(tuple(e) for e in existing_edges_arr)
    frac_pos = {i + 1: frac_pos_arr[i] for i in range(len(frac_pos_arr))}
    pocket_radii = {i + 1: float(pocket_radii_arr[i]) for i in range(len(pocket_radii_arr))}
    base_labels = sorted(frac_pos.keys())

    print(f"\nBase-cell nodes: {base_labels}")
    print(f"Cell: {cellpar}")
    print(f"Existing edges: {len(existing_set)}")

    # ── Step 2: enumerate all candidate offsets & build distance table ─
    offset_range = [-1, 0, 1]
    # Group existing edges by (i,j) pair
    pair_edges = {}
    for (i, j, off) in existing_set:
        pair_edges.setdefault((i, j), set()).add(off)

    print(f"\nNode pairs with edges: {len(pair_edges)}")

    # ── Step 3: distance-shell analysis per pair ──────────────────────
    all_missing = []

    for (i, j), existing_offsets in sorted(pair_edges.items()):
        fi, fj = frac_pos[i], frac_pos[j]
        candidates = []

        for di in offset_range:
            for dj in offset_range:
                for dk in offset_range:
                    offset = (di, dj, dk)
                    key = canonical_key(i, j, offset)
                    if key is None:
                        continue
                    dist, d_cart = cart_distance(fi, fj, np.array(offset), cell)
                    candidates.append((key, dist, d_cart))

        # Cluster into shells
        dist_list = [(key, dist) for key, dist, _ in candidates]
        shells = cluster_shell(dist_list, tol=0.15)

        for shell in shells:
            shell_keys = {k for k, _ in shell}
            shell_exists = shell_keys & existing_set
            shell_missing = shell_keys - existing_set

            # Only flag if the shell is *partially* populated
            if shell_exists and shell_missing:
                ref_dist = shell[0][1]
                for key, dist in shell:
                    if key in shell_missing:
                        all_missing.append({
                            "edge": key,
                            "dist": dist,
                            "shell_dist": ref_dist,
                            "shell_size": len(shell),
                            "shell_existing": len(shell_exists),
                            "shell_missing": len(shell_missing),
                            "node_i": key[0],
                            "node_j": key[1],
                            "offset": key[2],
                        })

    # ── Step 4: report ────────────────────────────────────────────────
    print(f"\n{'=' * 72}")
    print(f"MISSING EDGE DETECTION RESULTS")
    print(f"{'=' * 72}")
    print(f"Total missing edges found: {len(all_missing)}")

    # Group by node pair
    by_pair = {}
    for m in all_missing:
        pair = (m["node_i"], m["node_j"])
        by_pair.setdefault(pair, []).append(m)

    for pair, missing_list in sorted(by_pair.items()):
        i, j = pair
        existing = pair_edges.get(pair, set())
        print(f"\n--- Pair ({i}, {j}) ---")
        print(f"  Pocket radii: r_{i}={pocket_radii[i]:.2f}, r_{j}={pocket_radii[j]:.2f}")
        print(f"  Existing offsets ({len(existing)}):")
        for off in sorted(existing):
            dist, _ = cart_distance(frac_pos[i], frac_pos[j], np.array(off), cell)
            print(f"    ({i}, {j}, {off})  |d|={dist:.2f} Å")
        print(f"  Missing offsets ({len(missing_list)}):")
        for m in sorted(missing_list, key=lambda x: x["dist"]):
            print(f"    ({i}, {j}, {m['offset']})  |d|={m['dist']:.2f} Å  "
                  f"[shell: {m['shell_existing']} exist + {m['shell_missing']} missing, "
                  f"ref={m['shell_dist']:.2f} Å]")

    # ── Step 5: save results ──────────────────────────────────────────
    result_path = Path("/tmp/nu1000_missing_edges_results.npz")
    np.savez(
        result_path,
        missing_edges=np.array([m["edge"] for m in all_missing], dtype=object),
        missing_dists=np.array([m["dist"] for m in all_missing]),
        missing_details=np.array(all_missing, dtype=object),
    )
    print(f"\nResults saved to {result_path}")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
