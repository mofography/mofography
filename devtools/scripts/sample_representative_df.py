"""
Stratified (representative) sample from a MOF dataframe so the sample has the same
rough distribution as the full dataset over: surface_area_m2g, void_fraction, lcd, pld.

Uses mofography.data.sample_representative_mofdb and plot_representative_sample_subplots.
Loads a feather dataframe from ~/Datasets (default: coredf.feather), saves the sample,
and optionally opens a Plotly figure in the browser.

Usage (from repo root):
  python devtools/scripts/sample_representative_df.py [--path PATH] [--n N] [--bins B] [--seed S]
  Default: path=~/Datasets/coredf.feather, n=2000, bins=4, seed=0
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd

import mofography as mgr


def load_core_df(path: str | Path):
    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Dataset not found: {path}")
    return pd.read_feather(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
    parser.add_argument(
        "--path",
        type=str,
        default=os.path.expanduser("~/Datasets/coredf.feather"),
        help="Path to feather dataframe",
    )
    parser.add_argument(
        "--n",
        type=int,
        default=2000,
        help="Target sample size",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=4,
        help="Number of bins per variable for stratification",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--out",
        type=str,
        default=None,
        help="Output path for sample (default: same dir as path, name with _sample_n<N>.feather)",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not open figure in browser (only save HTML if --save-html given)",
    )
    parser.add_argument(
        "--save-html",
        type=str,
        default=None,
        help="Save figure to this HTML path",
    )
    args = parser.parse_args()

    path = Path(args.path).expanduser().resolve()
    df = load_core_df(path)
    print(f"Loaded {path.name}: {len(df)} rows")

    sample = mgr.sample_representative_mofdb(
        df,
        n=args.n,
        strat_cols=mgr.MOFDB_STRAT_COLS,
        n_bins=args.bins,
        random_state=args.seed,
    )
    print(f"Sampled {len(sample)} materials (stratified, n_bins={args.bins})")

    out_path = args.out
    if out_path is None:
        out_path = path.parent / f"{path.stem}_sample_n{len(sample)}.feather"
    else:
        out_path = Path(out_path).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sample.to_feather(out_path)
    print(f"Saved sample to {out_path}")

    title = f"Representative Subset of CoRE MOFs 2019 (N={len(sample)})"
    fig = mgr.plot_representative_sample_subplots(
        df, sample, n_sample=len(sample), title=title
    )
    if args.save_html:
        fig.write_html(args.save_html)
        print(f"Saved figure to {args.save_html}")
    if not args.no_browser:
        fig.show()
    return sample, fig


if __name__ == "__main__":
    main()
