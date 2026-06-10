import sys
sys.path.insert(0, "/Users/arungopalan/repos/mofography")
import mofography as mgr
from pathlib import Path

base = Path("/Users/arungopalan/repos/mofography/outputs")
out_html = base / "github_raspa_test" / "flexible_irmof1_benzene_sidecar_test4.html"

out_path, payload = mgr.make_sparse_multicomponent_movie(
    component_movie_paths={
        "benzene": base / "raspa3_md_benzene_eval" / "movies" / "movie.s0.pdb",
    },
    framework_movie_path=base / "github_raspa_test" / "Movie_IRMOF1_frameworks.pdb",
    cif_path=base / "github_raspa_test" / "IRMOF-1.cif",
    out_html=out_html,
    max_frames=4,
    movie_format="raspa_pdb",
    spacing=0.25,
    isovalue=1.0,
    surf_step_size=4,
    sidecar=True,
)

print(f"Done! File: {out_path}")
print(f"Size: {out_path.stat().st_size / (1024*1024):.1f} MB")

sidecar = out_path.with_suffix(".frames.jsonl")
if sidecar.exists():
    print(f"Sidecar: {sidecar}")
    print(f"Sidecar size: {sidecar.stat().st_size / (1024*1024):.1f} MB")
    with open(sidecar) as f:
        lines = sum(1 for _ in f)
    print(f"Sidecar frames: {lines}")
else:
    print("No sidecar file found")
