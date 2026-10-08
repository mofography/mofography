# Releases and citing mofography

## Policy (for collaborators)

Cite a **released tag** `YYYY.MM.DD.N` from the [mofography releases](https://github.com/mofography/mofography/releases) page, not an arbitrary checkout of `master`. For reproducibility, record the **tag name** and the **full git commit hash** shown on that release.

## Version format

- Tags look like `2026.10.09.1` (year, month, day, same-day increment).
- Second release on the same calendar day: `2026.10.09.2`, and so on.
- Tags are **annotated** git tags on the commit you want others to reproduce.
- [versioneer](https://github.com/python-versioneer/python-versioneer) (`setup.cfg`, PEP 440) derives `mofography.__version__` from the nearest tag.

## Release checklist

1. Merge or commit everything that should be in the release; working tree clean (`git status` shows no modifications you care about).
2. Run tests you rely on for the paper or milestone.
3. Update `CITATION.cff`: set `version` to the new tag string (e.g. `2026.10.09.1`).
4. Commit the `CITATION.cff` bump (if not already on the release commit).
5. Create an annotated tag on that commit:

   ```bash
   git tag -a 2026.10.09.1 -m "mofography 2026.10.09.1 — <short reason>"
   ```

6. Push the commit and tag to the canonical remote:

   ```bash
   git push org master
   git push org 2026.10.09.1
   ```

7. Open a **GitHub Release** from the tag; paste a short note (paper, breaking changes, install hint).
8. Verify from a clean checkout at the tag:

   ```bash
   git checkout 2026.10.09.1
   python -c "import mofography as m; print(m.__version__, m.__git_revision__)"
   ```

   Expect `__version__` to match the tag (not `0+untagged…` or `.dirty`).

## What to put in a paper or SI

- Software name: **mofography**
- Version: tag string, e.g. `2026.10.09.1`
- URL: `https://github.com/mofography/mofography`
- Commit: full hash from the release (or `mofography.__git_revision__` at that tag)
- Optional: copy the metadata block from `CITATION.cff` after updating `version`

## Interactive HTML demos

The repo does **not** ship large Three.js HTML snapshots under `examples/viz_outputs/` (they go stale quickly). Use the mofography web app (local or public demo on Hugging Face) to build and download current HTML, or regenerate with the `examples/gen_*_threejs.py` scripts when you need a file on disk.

## Zenodo (later)

When the repository is public, enable [GitHub–Zenodo integration](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content) so each GitHub Release can receive a DOI. Until then, the git tag plus commit hash is the authoritative citeable snapshot.
