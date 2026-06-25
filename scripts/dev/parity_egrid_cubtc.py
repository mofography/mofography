"""Parity plots: Python egrid_from_atoms2 vs Rust scalar vs Rust tricubic on CuBTC.

Energies are in Kelvin (LJ epsilon / k_B). Only voxels with mutual |E| < 100 K are plotted.

Grid alignment at the same array index (i, j, k):
  - Rust: fractional node  frac_r = -0.5 + i/(N-1)  with MIC on (scaled - 0.5) atoms
  - Python (freud/AABB): ASE fractional frac_a = i/(N-1)  -> Cartesian via the unit cell
  These sample the same physical point because frac_a = frac_r + 0.5.
"""
import sys
import time
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from ase.io import read

_here = Path(__file__).resolve()
if (_here.parents[2] / "mofography" / "__init__.py").exists():
    ROOT = _here.parents[2]
else:
    ROOT = _here.parents[2] / "mofography"
sys.path.insert(0, str(ROOT))

import mofography as mgr

CIF = ROOT / "mofography/examples/TestCIFs/P1-Cu-BTC.cif"
SPACING = 0.5
CUTOFF = 12.0
PROBE = "He"
E_MAX_K = 100.0
OUT = Path("/Users/arungopalan/repos/mgr_sn/mgr_org_scripts/parity_outputs_cubtc_egrid")
OUT_MOFO = ROOT / "parity_outputs"
for _out in (OUT, OUT_MOFO):
    _out.mkdir(parents=True, exist_ok=True)

MARKER = dict(
    size=10,
    color="rgba(70, 130, 180, 0.65)",
    line=dict(width=1.4, color="rgb(25, 55, 90)"),
)


def _masked_pair(ref, comp, e_max=E_MAX_K):
    mask = (
        np.isfinite(ref)
        & np.isfinite(comp)
        & (ref < e_max)
        & (comp < e_max)
        & (ref > -e_max)
        & (comp > -e_max)
    )
    return ref[mask].ravel(), comp[mask].ravel(), mask.sum()


def _subsample(x, y, max_points=25_000, seed=0):
    if x.size <= max_points:
        return x, y
    rng = np.random.default_rng(seed)
    idx = rng.choice(x.size, max_points, replace=False)
    return x[idx], y[idx]


def _parity_trace(xp, yp, ref_label, comp_label, name="grid points"):
    abs_dev = np.abs(yp - xp)
    return go.Scatter(
        x=xp,
        y=yp,
        mode="markers",
        name=name,
        customdata=np.column_stack([yp, abs_dev]),
        hovertemplate=(
            f"<b>{ref_label}</b>: %{{x:.4f}} K<br>"
            f"<b>{comp_label}</b>: %{{customdata[0]:.4f}} K<br>"
            "<b>|Δ|</b>: %{customdata[1]:.4f} K"
            "<extra></extra>"
        ),
        marker=MARKER,
    )


def _identity_trace(ax_lo, ax_hi, name="y = x", showlegend=False):
    return go.Scatter(
        x=[ax_lo, ax_hi],
        y=[ax_lo, ax_hi],
        mode="lines",
        line=dict(color="crimson", dash="dash", width=2),
        name=name,
        showlegend=showlegend,
    )


def _stats(ref, comp):
    diff = comp - ref
    return {
        "rmse": float(np.sqrt(np.mean(diff**2))),
        "mad": float(np.mean(np.abs(diff))),
        "max_abs": float(np.max(np.abs(diff))),
        "n": int(ref.size),
        "r2": float(np.corrcoef(ref, comp)[0, 1] ** 2),
    }


def make_parity(ref, comp, ref_label, comp_label, title, outpath, e_max=E_MAX_K):
    """Scatter parity on points with both energies < e_max K."""
    x, y, _ = _masked_pair(ref, comp, e_max=e_max)
    if x.size == 0:
        raise RuntimeError(f"No points below {e_max} K for {title}")

    st = _stats(x, y)
    xp, yp = _subsample(x, y)

    lo = float(min(x.min(), y.min()))
    hi = float(max(x.max(), y.max()))
    pad = 0.02 * (hi - lo + 1e-9)
    ax_lo, ax_hi = lo - pad, hi + pad

    fig = go.Figure()
    fig.add_trace(_parity_trace(xp, yp, ref_label, comp_label))
    fig.add_trace(_identity_trace(ax_lo, ax_hi, showlegend=True))
    fig.update_layout(
        title=(
            f"{title}<br>"
            f"<sup>n = {st['n']:,} points with |E| &lt; {e_max:.0f} K · "
            f"R² = {st['r2']:.8f} · RMSE = {st['rmse']:.4f} K · "
            f"MAD = {st['mad']:.4f} K · max|Δ| = {st['max_abs']:.4f} K</sup>"
        ),
        xaxis_title=f"{ref_label} energy (K)",
        yaxis_title=f"{comp_label} energy (K)",
        xaxis=dict(range=[ax_lo, ax_hi], scaleanchor="y", scaleratio=1),
        yaxis=dict(range=[ax_lo, ax_hi]),
        width=620,
        height=620,
        template="plotly_white",
        showlegend=True,
    )
    fig.write_html(str(outpath))
    print(f"  → {outpath}")
    return st


def make_combined_parity(
    egrid_py,
    egrid_rust,
    egrid_tri,
    outpaths,
    e_max=E_MAX_K,
):
    """One HTML figure with three parity subplots: Rust/Python, Tricubic/Python, Rust/Tricubic."""
    panels = [
        (
            egrid_py,
            egrid_rust,
            "Python egrid_from_atoms2",
            "Rust egrid_from_atoms",
            "Rust scalar vs Python",
        ),
        (
            egrid_py,
            egrid_tri,
            "Python egrid_from_atoms2",
            "Rust tricubic (ch0)",
            "Rust tricubic vs Python",
        ),
        (
            egrid_rust,
            egrid_tri,
            "Rust egrid_from_atoms",
            "Rust tricubic (ch0)",
            "Rust tricubic vs Rust scalar",
        ),
    ]

    stats = []
    fig = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=[p[4] for p in panels],
        horizontal_spacing=0.08,
    )

    global_lo, global_hi = np.inf, -np.inf
    for ref, comp, ref_label, comp_label, _title in panels:
        x, y, _ = _masked_pair(ref, comp, e_max=e_max)
        stats.append(_stats(x, y))
        global_lo = min(global_lo, float(x.min()), float(y.min()))
        global_hi = max(global_hi, float(x.max()), float(y.max()))

    pad = 0.02 * (global_hi - global_lo + 1e-9)
    ax_lo, ax_hi = global_lo - pad, global_hi + pad

    for col, (panel, st) in enumerate(zip(panels, stats), start=1):
        ref, comp, ref_label, comp_label, title = panel
        x, y, _ = _masked_pair(ref, comp, e_max=e_max)
        xp, yp = _subsample(x, y, seed=col)

        fig.add_trace(
            _parity_trace(xp, yp, ref_label, comp_label, name=title),
            row=1,
            col=col,
        )
        fig.add_trace(
            _identity_trace(ax_lo, ax_hi, showlegend=(col == 1)),
            row=1,
            col=col,
        )
        fig.update_xaxes(
            title_text=f"{ref_label} (K)",
            range=[ax_lo, ax_hi],
            scaleanchor=f"y{col if col > 1 else ''}",
            scaleratio=1,
            row=1,
            col=col,
        )
        fig.update_yaxes(
            title_text=f"{comp_label} (K)",
            range=[ax_lo, ax_hi],
            row=1,
            col=col,
        )
        fig.layout.annotations[col - 1].update(
            text=(
                f"<b>{title}</b><br>"
                f"<sup>n={st['n']:,} · R²={st['r2']:.6f} · "
                f"RMSE={st['rmse']:.4f} K · MAD={st['mad']:.4f} K</sup>"
            )
        )

    n_mutual = int(
        (
            np.isfinite(egrid_py)
            & np.isfinite(egrid_rust)
            & np.isfinite(egrid_tri)
            & (np.abs(egrid_py) < e_max)
            & (np.abs(egrid_rust) < e_max)
            & (np.abs(egrid_tri) < e_max)
        ).sum()
    )

    fig.update_layout(
        title=dict(
            text=(
                f"<b>CuBTC egrid parity (He probe, spacing {SPACING} Å, cutoff {CUTOFF} Å)</b><br>"
                f"<sup>Mutual |E| &lt; {e_max:.0f} K · {n_mutual:,} voxels · energies in K (ε/k<sub>B</sub>)</sup>"
            ),
            x=0.5,
            xanchor="center",
        ),
        width=1900,
        height=640,
        template="plotly_white",
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=-0.18, x=0.5, xanchor="center"),
    )

    for outpath in outpaths:
        fig.write_html(str(outpath))
        print(f"  → {outpath}")
    return stats


def make_diff_hist(ref, comp, ref_label, comp_label, title, outpath, e_max=E_MAX_K):
    mask = (
        np.isfinite(ref)
        & np.isfinite(comp)
        & (ref < e_max)
        & (comp < e_max)
        & (ref > -e_max)
        & (comp > -e_max)
    )
    diff = (comp[mask].ravel() - ref[mask].ravel())
    st = _stats(ref[mask].ravel(), comp[mask].ravel())
    lo, hi = np.percentile(diff, [0.5, 99.5])

    fig = go.Figure()
    fig.add_trace(
        go.Histogram(x=diff, nbinsx=120, marker=dict(color="steelblue"), name="ΔE")
    )
    fig.update_layout(
        title=(
            f"{title}<br>"
            f"<sup>ΔE = ({comp_label}) − ({ref_label}) · "
            f"mean = {diff.mean():.4f} K · std = {diff.std():.4f} K · "
            f"RMSE = {st['rmse']:.4f} K · MAD = {st['mad']:.4f} K</sup>"
        ),
        xaxis_title="Energy difference (K)",
        yaxis_title="Count",
        xaxis_range=[float(lo), float(hi)],
        width=700,
        height=420,
        template="plotly_white",
    )
    fig.write_html(str(outpath))
    print(f"  → {outpath}")


def main():
    print(f"Loading {CIF}...")
    atoms = read(str(CIF))
    cellpar = atoms.cell.cellpar()
    print(
        f"  {len(atoms)} atoms, cell "
        f"{cellpar[0]:.1f}×{cellpar[1]:.1f}×{cellpar[2]:.1f} Å"
    )

    print("\n[1/3] Python egrid_from_atoms2...")
    t0 = time.time()
    egrid_py = np.asarray(
        mgr.egrid_from_atoms2(
            atoms, spacing=SPACING, cutoff=CUTOFF, probe_symbol=PROBE
        ).compute(),
        dtype=np.float64,
    )
    t_py = time.time() - t0

    print("[2/3] Rust egrid_from_atoms_rust...")
    t0 = time.time()
    egrid_rust = np.asarray(
        mgr.egrid_from_atoms_rust(
            atoms, spacing=SPACING, cutoff=CUTOFF, probe_symbol=PROBE
        ),
        dtype=np.float64,
    )
    t_rust = time.time() - t0

    print("[3/3] Rust egrid_tricubic_from_atoms_rust (channel 0)...")
    t0 = time.time()
    egrid_tri = np.asarray(
        mgr.egrid_tricubic_from_atoms_rust(
            atoms, spacing=SPACING, cutoff=CUTOFF, probe_symbol=PROBE
        )[0],
        dtype=np.float64,
    )
    t_tri = time.time() - t0

    print(
        f"\nShapes: py={egrid_py.shape}, rust={egrid_rust.shape}, tri={egrid_tri.shape}"
    )
    for name, arr in [
        ("Python", egrid_py),
        ("Rust scalar", egrid_rust),
        ("Rust tricubic", egrid_tri),
    ]:
        finite = np.isfinite(arr)
        print(
            f"  {name}: min={arr[finite].min():.2f} K, max={arr[finite].max():.2e} K"
        )

    mutual = (
        np.isfinite(egrid_py)
        & np.isfinite(egrid_rust)
        & np.isfinite(egrid_tri)
        & (np.abs(egrid_py) < E_MAX_K)
        & (np.abs(egrid_rust) < E_MAX_K)
        & (np.abs(egrid_tri) < E_MAX_K)
    )
    print(f"\nMutual |E| < {E_MAX_K:.0f} K: {mutual.sum():,} / {egrid_py.size:,} voxels")

    print("\n--- Parity plots (|E| < 100 K) ---")
    s_pr = make_parity(
        egrid_py,
        egrid_rust,
        "Python egrid_from_atoms2",
        "Rust egrid_from_atoms",
        "CuBTC · Rust scalar vs Python",
        OUT / "parity_rust_vs_python.html",
    )
    make_diff_hist(
        egrid_py,
        egrid_rust,
        "Python egrid_from_atoms2",
        "Rust egrid_from_atoms",
        "CuBTC · Rust − Python difference",
        OUT / "diff_rust_vs_python.html",
    )

    s_pt = make_parity(
        egrid_py,
        egrid_tri,
        "Python egrid_from_atoms2",
        "Rust tricubic (ch0)",
        "CuBTC · Rust tricubic vs Python",
        OUT / "parity_tricubic_vs_python.html",
    )
    make_diff_hist(
        egrid_py,
        egrid_tri,
        "Python egrid_from_atoms2",
        "Rust tricubic (ch0)",
        "CuBTC · Tricubic − Python difference",
        OUT / "diff_tricubic_vs_python.html",
    )

    s_rt = make_parity(
        egrid_rust,
        egrid_tri,
        "Rust egrid_from_atoms",
        "Rust tricubic (ch0)",
        "CuBTC · Rust tricubic vs Rust scalar",
        OUT / "parity_rust_vs_tricubic.html",
    )

    print("\n--- Combined parity figure ---")
    combined_stats = make_combined_parity(
        egrid_py,
        egrid_rust,
        egrid_tri,
        outpaths=[
            OUT / "parity_cubtc_egrid_combined.html",
            OUT_MOFO / "parity_cubtc_egrid_combined.html",
        ],
    )
    s_pr, s_pt, s_rt = combined_stats[0], combined_stats[1], combined_stats[2]

    print(f"\n{'=' * 60}")
    print(f"Summary (mutual |E| < {E_MAX_K:.0f} K, energies in K):")
    for label, st in [
        ("Rust scalar vs Python", s_pr),
        ("Rust tricubic vs Python", s_pt),
        ("Rust scalar vs tricubic", s_rt),
    ]:
        print(
            f"  {label}: R²={st['r2']:.8f}, "
            f"RMSE={st['rmse']:.6f} K, MAD={st['mad']:.6f} K, "
            f"max|Δ|={st['max_abs']:.6f} K"
        )
    print(f"  Timing: Python={t_py:.1f}s, Rust={t_rust:.1f}s, Tricubic={t_tri:.1f}s")
    print(f"\nOutputs: {OUT}/")


if __name__ == "__main__":
    main()
