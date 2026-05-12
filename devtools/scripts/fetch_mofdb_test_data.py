"""
Fetch a small MOFDB slice for testing (IZA or a small CoREMOF subset)
and save as compressed parquet under ~/Datasets/.

Usage (from repo root):
  python devtools/scripts/fetch_mofdb_test_data.py [--iza]
  Default: small CoREMOF 2019 slice (vf 0.7–0.725, lcd 12–14)
"""
import argparse
import asyncio
import os
from pathlib import Path

import mofography as mgr


def get_test_params(use_iza: bool):
    if use_iza:
        return dict(vf_min=0.2, vf_max=1.0, database="IZA")
    return dict(
        vf_min=0.7,
        vf_max=0.725,
        lcd_min=12,
        lcd_max=14,
        database="CoREMOF 2019",
    )


async def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split("\n")[0])
    parser.add_argument(
        "--iza",
        action="store_true",
        help="Use IZA database (small). Default is small CoREMOF 2019 subset.",
    )
    args = parser.parse_args()
    use_iza = args.iza

    params = get_test_params(use_iza=use_iza)
    df = await mgr.get_data_from_mofdb(params=params)

    out_dir = Path(os.path.expanduser("~/Datasets"))
    out_dir.mkdir(parents=True, exist_ok=True)
    if use_iza:
        out_path = out_dir / "iza_test.parquet"
    else:
        out_path = out_dir / "CoREMOF_2019_subset_test.parquet"
    df.to_parquet(out_path, compression="zstd")
    print(f"Saved {len(df)} rows to {out_path}")


if __name__ == "__main__":
    asyncio.run(main())
