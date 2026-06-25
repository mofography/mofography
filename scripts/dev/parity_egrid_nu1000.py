"""Parity: Rust scalar egrid vs tricubic *interpolator* on NU-1000 grid nodes.

Non-orthogonal cell (ASE cell matrix rows; Rust uses cellpar_to_cell(...).T).

At each grid index (i, j, k) the scalar energy is egrid_from_atoms_rust[i,j,k].
The tricubic value is TricubicField.interpolate(frac_ase) at the same physical
node — NOT channel 0 of the 8-grid directly.

Grid node ASE fractional coords (matches Rust kernel):
  frac_ase_a = (i / (N_a - 1)) * extent_a
with extent from repeat_to_meet_cutoff (1.0 when no replication).

Only points with mutual |E| < 100 K are plotted (Plotly HTML).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from ase.io import read

_here = Path(__file__).resolve()
ROOT = _here.parents[2] if (_here.parents[2] / "mofography" / "__init__.py").exists() else _here.parents[2] / "mofography"
sys.path.insert(0, str(ROOT))

import mofography as mgr
from mofography.rust_fun import (
    egrid_query_extent_from_atoms,
    egrid_tricubic_from_atoms_rust,
    tricubic_interpolator_from_grid,
)

CIF = ROOT / "examples/TestCIFs/NU-1000_S1.cif"
if not CIF.exists():
    CIF = ROOT / "mofography/examples/TestCIFs/NU-1000_S1.cif"
SPACING = 0.5
CUTOFF = 12.5
PROBE = "He"
E_MAX_K = 100.0
OUT = ROOT / "parity_outputs"
OUT_MGR = Path("/Users/arungopalan/repos/mgr_sn/mgr_org_scripts/parity_outputs_nu1000_egrid")
for _out in (OUT, OUT_MGR):
    _out.mkdir(parents=True, exist_ok=True)

MARKER = dict(
    size=11,
    color="rgba(46, 139, 87, 0.7)",
    line=dict(width=1.5, color="rgb(20, 70, 40)"),
)


def grid_node_frac_ase(nx: int, ny: int, nz: int, extent) -> np.ndarray:
    """(N, 3) ASE fractional coords for each grid voxel center."""
    ex, ey, ez = extent
    ix = np.arange(nx, dtype=np.float64) / max(nx - 1, 1)
    iy = np.arange(ny, dtype=np.float64) / max(ny - 1, 1)
    iz = np.arange(nz, dtype=np.float64) / max(nz - 1, 1)
    fx = ix * ex
    fy = iy * ey
    fz = iz * ez
    ii, jj, kk = np.meshgrid(fx, fy, fz, indexing="ij")
    return np.column_stack([ii.ravel(), jj.ravel(), kk.ravel()])


def _masked_pair(ref, comp, e_max=E_MAX_K):
    mask = (
        np.isfinite(ref)
        & np.isfinite(comp)
        & (np.abs(ref) < e_max)
        & (np.abs(comp) < e_max)
    )
    return ref[mask].ravel(), comp[mask].ravel(), int(mask.sum())


def _subsample(x, y, max_points=30_000, seed=0):
    if x.size <= max_points:
        return x, y
    rng = np.random.default_rng(seed)
    idx = rng.choice(x.size, max_points, replace=False)
    return x[idx], y[idx]


def _stats(ref, comp):
    diff = comp - ref
    return {
        "rmse": float(np.sqrt(np.mean(diff**2))),
        "mad": float(np.mean(np.abs(diff))),
        "max_abs": float(np.max(np.abs(diff))),
        "n": int(ref.size),
        "r2": float(np.corrcoef(ref, comp)[0, 1] ** 2) if ref.size > 1 else float("nan"),
    }


def make_parity(ref, comp, ref_label, comp_label, title, outpath, e_max=E_MAX_K):
    x, y, n = _masked_pair(ref, comp, e_max=e_max)
    if x.size == 0:
        raise RuntimeError(f"No mutual |E| < {e_max} K points for {title}")
    st = _stats(x, y)
    xp, yp = _subsample(x, y)
    abs_dev = np.abs(yp - xp)

    lo = float(min(x.min(), y.min()))
    hi = float(max(x.max(), y.max()))
    pad = 0.02 * (hi - lo + 1e-9)
    ax_lo, ax_hi = lo - pad, hi + pad

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=xp,
            y=yp,
            mode="markers",
            name="grid nodes",
            customdata=np.column_stack([yp, abs_dev]),
            hovertemplate=(
                f"<b>{ref_label}</b>: %{{x:.4f}} K<br>"
                f"<b>{comp_label}</b>: %{{customdata[0]:.4f}} K<br>"
                "<b>|Δ|</b>: %{customdata[1]:.4f} K"
                "<extra></extra>"
            ),
            marker=MARKER,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[ax_lo, ax_hi],
            y=[ax_lo, ax_hi],
            mode="lines",
            line=dict(color="crimson", dash="dash", width=2),
            name="y = x",
        )
    )
    fig.update_layout(
        title=(
            f"{title}<br>"
            f"<sup>n = {st['n']:,} nodes with |E| &lt; {e_max:.0f} K · "
            f"R² = {st['r2']:.8f} · RMSE = {st['rmse']:.4f} K · "
            f"MAD = {st['mad']:.4f} K · max|Δ| = {st['max_abs']:.4f} K</sup>"
        ),
        xaxis_title=f"{ref_label} energy (K)",
        yaxis_title=f"{comp_label} energy (K)",
        xaxis=dict(range=[ax_lo, ax_hi], scaleanchor="y", scaleratio=1),
        yaxis=dict(range=[ax_lo, ax_hi]),
        width=680,
        height=680,
        template="plotly_white",
    )
    fig.write_html(str(outpath))
    print(f"  → {outpath}")
    return st


def main():
    print(f"Loading {CIF}...")
    atoms = read(str(CIF))
    cp = atoms.cell.cellpar()
    print(
        f"  {len(atoms)} atoms, cell "
        f"{cp[0]:.2f}×{cp[1]:.2f}×{cp[2]:.2f} Å, "
        f"α={cp[3]:.1f} β={cp[4]:.1f} γ={cp[5]:.1f}°"
    )

    extent = egrid_query_extent_from_atoms(atoms, cutoff=CUTOFF)
    print(f"  egrid query extent (ASE tile): {extent}")

    print("\n[1/2] Rust scalar egrid_from_atoms_rust...")
    t0 = time.time()
    egrid = np.asarray(
        mgr.egrid_from_atoms_rust(
            atoms, spacing=SPACING, cutoff=CUTOFF, probe_symbol=PROBE
        ),
        dtype=np.float64,
    )
    t_scalar = time.time() - t0
    nx, ny, nz = egrid.shape
    print(f"  shape {egrid.shape}, build {t_scalar:.1f}s")

    print("[2/2] Rust tricubic grid + interpolator...")
    t0 = time.time()
    egrid8 = np.asarray(
        egrid_tricubic_from_atoms_rust(
            atoms, spacing=SPACING, cutoff=CUTOFF, probe_symbol=PROBE
        ),
        dtype=np.float64,
    )
    field = tricubic_interpolator_from_grid(egrid8, *extent)
    t_tri_build = time.time() - t0

    frac_nodes = grid_node_frac_ase(nx, ny, nz, extent)
    print(f"  interpolating at {frac_nodes.shape[0]:,} nodes...")
    t0 = time.time()
    e_tri_interp = field.interpolate_batch(frac_nodes).reshape(nx, ny, nz)
    t_interp = time.time() - t0

    e_ch0 = egrid8[0]
    finite = np.isfinite(egrid)
    print(
        f"  scalar: min={egrid[finite].min():.2f} K, max={egrid[finite].max():.2e} K"
    )
    print(
        f"  tricubic interp: min={e_tri_interp[finite].min():.2f} K, "
        f"max={e_tri_interp[finite].max():.2e} K"
    )

    mutual_interp = (
        np.isfinite(egrid)
        & np.isfinite(e_tri_interp)
        & (np.abs(egrid) < E_MAX_K)
        & (np.abs(e_tri_interp) < E_MAX_K)
    ).sum()
    print(f"\nMutual |E| < {E_MAX_K:.0f} K (scalar vs interp): {mutual_interp:,}")

    print("\n--- Parity plots ---")
    st_main = make_parity(
        egrid,
        e_tri_interp,
        "Rust scalar egrid",
        "Tricubic interpolator",
        f"NU-1000 · scalar vs tricubic interp (He, {SPACING} Å, cutoff {CUTOFF} Å)",
        OUT / "parity_nu1000_scalar_vs_tricubic_interp.html",
    )
    make_parity(
        egrid,
        e_tri_interp,
        "Rust scalar egrid",
        "Tricubic interpolator",
        f"NU-1000 · scalar vs tricubic interp",
        OUT_MGR / "parity_nu1000_scalar_vs_tricubic_interp.html",
    )

    st_ch0 = make_parity(
        egrid,
        e_ch0,
        "Rust scalar egrid",
        "Tricubic grid ch0",
        f"NU-1000 · scalar vs tricubic ch0 (build check)",
        OUT / "parity_nu1000_scalar_vs_tricubic_ch0.html",
    )

    print(f"\n{'=' * 60}")
    print(f"Scalar vs tricubic interpolator: R²={st_main['r2']:.8f}, "
          f"RMSE={st_main['rmse']:.6f} K, MAD={st_main['mad']:.6f} K")
    print(f"Scalar vs tricubic ch0:         R²={st_ch0['r2']:.8f}, "
          f"RMSE={st_ch0['rmse']:.6f} K, MAD={st_ch0['mad']:.6f} K")
    print(f"Timing: scalar={t_scalar:.1f}s, tri build={t_tri_build:.1f}s, "
          f"interp batch={t_interp:.1f}s")
    print(f"\nOpen: {OUT / 'parity_nu1000_scalar_vs_tricubic_interp.html'}")


if __name__ == "__main__":
    main()
