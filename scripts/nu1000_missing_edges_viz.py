#!/usr/bin/env python
"""
NU-1000: detect missing close_or_inside edges, append to connections,
and generate enhanced HTML visualization with regions + regions_pbc + maxima.
"""

import sys
from pathlib import Path

import numpy as np

MOFOGRAPHY_ROOT = Path("/Users/arungopalan/repos/mofography")
sys.path.insert(0, str(MOFOGRAPHY_ROOT))
sys.path.insert(0, str(Path.home() / "repos/tipmof"))

from mofography.cif_io import read_cif_to_ase
import mofography as mgr
from tipmof.tipmof import collect_periodic_edges

CIF = str(MOFOGRAPHY_ROOT / "examples/TestCIFs/NU-1000.cif")
OUTPUT = "/tmp/nu1000_poregraph_with_missing.html"

# ── helpers ──

def canonical_key(i, j, offset):
    offset = tuple(int(x) for x in offset)
    if i < j: return (i, j, offset)
    elif i > j: return (j, i, tuple(-x for x in offset))
    else:
        if offset == (0,0,0): return None
        neg = tuple(-x for x in offset)
        return (i, i, offset) if offset < neg else (i, i, neg)


def cart_dist(fi, fj, off, cell):
    d = fj + np.array(off, float) - fi
    return np.linalg.norm(d @ cell.T)


def cluster_shell(items, tol=0.15):
    if not items: return []
    s = sorted(items, key=lambda x: x[1])
    shells, cur, ref = [], [s[0]], s[0][1]
    for k, d in s[1:]:
        if ref > 0 and abs(d - ref)/ref < tol:
            cur.append((k, d))
        else:
            shells.append(cur); cur, ref = [(k, d)], d
    shells.append(cur)
    return shells


def edge_to_ci_endpoints(base_i, base_j, offset, ci_lookup):
    """Return (ci_a, ci_b) if both endpoints have close_or_inside regions."""
    # canonical form i<j: i at origin, j at +offset
    ep_i = ci_lookup.get((base_i, (0,0,0)))
    ep_j = ci_lookup.get((base_j, offset))
    if ep_i is not None and ep_j is not None:
        return ep_i, ep_j
    # try reverse
    ep_i2 = ci_lookup.get((base_i, tuple(-x for x in offset)))
    ep_j2 = ci_lookup.get((base_j, (0,0,0)))
    if ep_i2 is not None and ep_j2 is not None:
        return ep_i2, ep_j2
    return None, None


# ── main pipeline ──

print("Step 1: Running segmentation …")
atoms = read_cif_to_ase(CIF)
atoms = mgr.update_with_covalent_bonds(atoms)
mgr.compute_aromatic_cycles(atoms)

dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=0.5)
seg = mgr.do_periodic_segmentation_rust(
    dgrid, atoms, spacing=0.5, mask_thickness=0.05, h=0.2,
    threshold_abs=1.5, min_distance=1.5, extra_padding=5,
    compute_connections=True, return_dgrid=True)

from ase.geometry import cellpar_to_cell
cell = cellpar_to_cell(atoms.cell.cellpar())
maxima_info = seg["maxima_info"]
grid_shape = dgrid.shape

# ── Build ci_region lookup ──

ci_maxima = np.array(seg["close_or_inside_maxima"])
ci_regions = np.array(seg["close_or_inside_regions"])
ci_lookup = {}  # (base_label, offset) → ci_region_id

for ci_pos, ci_rid in zip(ci_maxima, ci_regions):
    ci_pos = np.array(ci_pos, int)
    for row in maxima_info:
        base_vox = row[1:4].astype(int)
        base_label = int(row[0])
        diff = ci_pos - base_vox
        ok = True
        for ax in range(3):
            if grid_shape[ax] == 0: ok = False; break
            oa = diff[ax] / grid_shape[ax]
            if abs(oa - round(oa)) > 0.1: ok = False; break
        if ok:
            offset = tuple(int(round(diff[ax] / grid_shape[ax])) for ax in range(3))
            ci_lookup[(base_label, offset)] = int(ci_rid)
            break

print(f"  close_or_inside regions: {len(ci_lookup)}")
print(f"  ci_lookup: { {f'base{k[0]}@{k[1]}': v for k,v in sorted(ci_lookup.items())} }")

# ── Build pore graph and find existing edges ──

edge_map, _ = collect_periodic_edges(seg["connections"], maxima_info)
existing = set(edge_map.keys())

# Frac positions
is_base = np.all(maxima_info[:, 4:7] == 0, axis=1)
base_labels = sorted(int(maxima_info[i, 0]) for i in np.where(is_base)[0])
frac_pos = {}
for node in base_labels:
    row = maxima_info[maxima_info[:, 0] == node][0]
    frac_pos[node] = np.round(row[1:4] / np.array(grid_shape) + row[4:7], 3)

# ── Detect missing edges ──

pair_edges = {}
for (i, j, off) in existing:
    pair_edges.setdefault((i, j), set()).add(off)

missing_ci = []  # missing edges between ci regions
for (i, j), ex_offs in sorted(pair_edges.items()):
    fi, fj = frac_pos[i], frac_pos[j]
    ex_dists = [cart_dist(fi, fj, np.array(o), cell) for o in ex_offs]
    cutoff = max(ex_dists) * 1.05
    cands = []
    for di in [-1,0,1]:
        for dj in [-1,0,1]:
            for dk in [-1,0,1]:
                key = canonical_key(i, j, (di,dj,dk))
                if key is None: continue
                d = cart_dist(fi, fj, np.array(key[2]), cell)
                if d <= cutoff:
                    cands.append((key, d))
    shells = cluster_shell(cands, tol=0.15)
    for shell in shells:
        skeys = {k for k,_ in shell}
        sexist = skeys & existing
        smiss = skeys - existing
        if sexist and smiss:
            for k, d in shell:
                if k in smiss:
                    ci_a, ci_b = edge_to_ci_endpoints(k[0], k[1], k[2], ci_lookup)
                    if ci_a is not None and ci_b is not None:
                        missing_ci.append({
                            "edge": k, "dist": d,
                            "ci_a": ci_a, "ci_b": ci_b,
                            "shell_ex": len(sexist), "shell_miss": len(smiss),
                        })

print(f"\n  Missing edges between ci regions: {len(missing_ci)}")
for m in sorted(missing_ci, key=lambda x: x["dist"]):
    k = m["edge"]
    print(f"    ({k[0]:>2},{k[1]:>2},{k[2]})  |d|={m['dist']:.2f}Å  "
          f"ci_regions ({m['ci_a']},{m['ci_b']})  "
          f"[{m['shell_ex']}ex+{m['shell_miss']}mi in shell]")

# ── Append missing connections to seg dict ──

# Connection format: [labels_array, is_valid_bool, type_str, boundary_pts, window_radii]
# For missing edges, we synthesize a placeholder connection with the two ci_region labels
import copy
seg_aug = copy.deepcopy(seg)
conns_aug = list(seg_aug["connections"])

for m in missing_ci:
    ci_a, ci_b = m["ci_a"], m["ci_b"]
    # Create a synthetic connection [labels, valid, type, boundary_pts, window_radii]
    midpt_frac = (frac_pos[m["edge"][0]] + frac_pos[m["edge"][1]] + np.array(m["edge"][2], float)) / 2
    midpt_frac = midpt_frac % 1.0
    midpt_cart = midpt_frac @ cell.T
    # Estimate window radius as average of pocket radii
    # (use a reasonable placeholder since we don't have the real one)
    wr = 2.0  # placeholder
    new_conn = [
        np.array([ci_a, ci_b]),     # labels
        True,                         # is_valid
        "synthetic_missing",          # type (tagged for identification)
        np.array([midpt_cart]),       # boundary point at midpoint
        np.array([wr]),               # placeholder window radius
    ]
    conns_aug.append(np.array(new_conn, dtype=object))

seg_aug["connections"] = conns_aug
print(f"\n  Original connections: {len(seg['connections'])}")
print(f"  Augmented connections: {len(conns_aug)}")

# ── Generate HTML with surfaces + maxima ──

print(f"\nStep 2: Generating HTML visualization …")
html = mgr.show_mof_threejs(
    seg_aug, atoms,
    plot_framework_representation="scatter",
    plot_surface=True,
    plot_surface_step_size=2,
    plot_connections=True,
    plot_maxima=True,
    plot_windows=False,
    plot_cell=True,
    output_html=OUTPUT)

print(f"\nSaved to {html}")
