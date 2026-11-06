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

## Example wording (collaborators)

Use the **tag and commit from the release you actually ran**, not `master`. Below uses release `2026.10.09.1` (commit `e608d33`) as a template—swap version/hash when citing a newer release.

### Methods or software paragraph (journal)

> Pore segmentation, graph construction, and visualization were performed with **mofography** version **2026.10.09.1** (Gopalan; https://github.com/mofography/mofography; release tag `2026.10.09.1`, commit `e608d3385211e2431c8ff5702db761c50011eb51`).

### Supporting Information (one line)

> mofography v2026.10.09.1 — https://github.com/mofography/mofography/releases/tag/2026.10.09.1

### Proposal (facilities / software section)

> **Open-source software:** We use **mofography** (v2026.10.09.1), a Python package for pore analysis in nanoporous materials, obtained from the canonical repository https://github.com/mofography/mofography at release tag `2026.10.09.1` (commit `e608d33`). Interactive 3D outputs are generated with the project web application; archived analysis uses the tagged source code only.

### LaTeX (manual `\bibitem` or supplementary refs)

```latex
\bibitem{gopalan2026mofography}
A.~Gopalan,
\textit{mofography} (version 2026.10.09.1),
\url{https://github.com/mofography/mofography},
release tag 2026.10.09.1, commit e608d3385211e2431c8ff5702db761c50011eb51,
2026.
```

After Zenodo is enabled, add `doi = {...}` to that entry and cite the **version DOI** in addition to (or instead of) the raw GitHub URL, per journal rules.

### Where `CITATION.cff` fits

Collaborators do not need to edit `CITATION.cff` to write a paper. It is the machine-readable record of the same metadata:

- At a release: https://github.com/mofography/mofography/blob/2026.10.09.1/CITATION.cff  
- GitHub shows a “Cite this repository” widget on the repo home when `CITATION.cff` is on the default branch.

To generate BibTeX from the file locally: [Citation File Format tools](https://citation-file-format.github.io/) or `cffconvert -i CITATION.cff -f bibtex` (if `cffconvert` is installed).

## Interactive HTML demos

The repo does **not** ship large Three.js HTML snapshots under `examples/viz_outputs/` (they go stale quickly). Use the mofography web app (local or public demo on Hugging Face) to build and download current HTML, or regenerate with the `examples/gen_*_threejs.py` scripts when you need a file on disk.

## Zenodo reminder (do this when the paper / repo goes public)

**Until `mofography/mofography` is public:** cite **GitHub release tag + commit** (e.g. `2026.10.09.1`, `e608d33`). That is enough for the manuscript and for collaborators.

**After arXiv (or when you intentionally open the org repo):** archive releases on Zenodo so journals can cite a **DOI**. Zenodo does not replace tagging; it snapshots each **GitHub Release** you already cut.

### Checklist (owner — do not skip)

1. [ ] Make **https://github.com/mofography/mofography** public (only when you are ready).
2. [ ] Log in to [zenodo.org](https://zenodo.org) with GitHub; grant access if prompted.
3. [ ] Zenodo → **Account → GitHub** → enable sync → turn **on** repository `mofography/mofography`.
4. [ ] For release `2026.10.09.1` (or the release tied to the paper): either create a new GitHub Release (already done) or use Zenodo’s “Create DOI” for that release if the integration offers a retry.
5. [ ] Open the new Zenodo deposit → review title, authors, license (BSD-3-Clause), version string → **Publish**.
6. [ ] Copy the **version DOI** (and note the **concept DOI** for “all versions of mofography”).
7. [ ] Add the version DOI to the paper SI / software note; optionally add `doi:` to `CITATION.cff` on a **future** release tag.
8. [ ] Tell collaborators: “Paper code = mofography **2026.10.09.1**; archival DOI = …”

### How others cite later (not necessarily the same DOI forever)

- **Your paper:** cite the version you used (`2026.10.09.1`) and, once available, that version’s **Zenodo DOI**.
- **Future users:** cite **their** release tag (and DOI for that tag). Zenodo’s **concept DOI** points at the software family; the **version DOI** is what reproducibility needs.

References: [GitHub–Zenodo](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content), [Zenodo GitHub guide](https://help.zenodo.org/docs/github/).
