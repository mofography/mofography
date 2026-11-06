# Building the documentation site

The site is a [Quarto](https://quarto.org) project in `docs/`. The API reference is generated with [quartodoc](https://machow.github.io/quartodoc/).

## Prerequisites

- [Quarto](https://quarto.org/docs/get-started/) 1.8 or newer
- The `mofography` Python environment, with the Rust extension built (see `docs/install.qmd`)
- Docs extras: `pip install quartodoc jupyter plotly pyyaml`

Tell Quarto which Python to use, so notebooks run against the environment that has `mofography`:

```bash
export QUARTO_PYTHON=$(which python)
```

## Build and view locally

```bash
cd docs
quarto render            # whole site -> docs/_site (takes several minutes)
quarto preview           # live preview in the browser
quarto render quick_start.qmd   # a single page
```

Open `docs/_site/index.html` in a browser. Pages that use `freeze: true` or `freeze: auto` re-use saved outputs from `docs/_freeze/` when present. Delete a page's folder in `docs/_freeze/` to force it to re-execute.

## Regenerating the API reference

Only needed after adding, removing, or re-documenting functions:

```bash
python scripts/dev/merge_function_status.py        # repo root: merge function status labels
python scripts/dev/split_deprecated_reference.py   # repo root: regroup reference by status
cd docs
quartodoc build                                     # writes docs/reference/*.qmd and _sidebar.yml
python ../scripts/dev/tag_reference.py              # adds status badges and the filter bar
quarto render
```

## Large interactive files

Interactive 3D HTML files (`docs/interactive/*.html`) are not stored in this repo. They live in the `mofography-gallery-assets` repository. Clone it and copy its HTML files into `docs/interactive/` before rendering if you want the gallery to show them.

## Sharing without GitHub Pages

After `quarto render`, zip `docs/_site` and send it. The site works offline: unzip and open `index.html`.

## Publishing on GitHub Pages

Live site: <https://mofography.github.io/mofography/>, served from the **`gh-pages` branch** (root `/`) of
`mofography/mofography`. There is **no CI**: `gh-pages` holds the rendered site and is updated by hand,
from a machine you control.

**Do not use `quarto publish gh-pages`.** It always pushes to the remote named `origin`, which is often a
personal fork rather than the org repo.

### Which pages can be re-rendered where

Tutorial pages execute code (data, the Rust extension, long runs). Their outputs are cached in
`docs/_freeze/` and `docs/.jupyter_cache/`, which are **gitignored and absent from a fresh clone**.

- **Text-only pages** (no executed code: `index.qmd`, `publications.qmd`, `install.qmd`, ...): safe to
  re-render in a fresh clone with the procedure below.
- **Pages with code** (tutorials): render only on a machine that has `docs/_freeze/` and the Python
  environment. Never delete `_freeze/`; it is the source of the stored outputs.

### Procedure: update some pages

```bash
# 1. Edit the source in a clone of master; commit and push it to the org repo.
git clone https://github.com/mofography/mofography.git && cd mofography
$EDITOR docs/publications.qmd
git add -A && git commit -m "Docs: ..." && git push origin master

# 2. Put the CURRENT published site into docs/_site, so sitemap.xml and search.json stay complete
#    and the large interactive/*.html files are kept.
rm -rf docs/_site
git clone --branch gh-pages https://github.com/mofography/mofography.git docs/_site

# 3. Render ONLY the changed pages (a single-file render also updates search.json and sitemap.xml).
cd docs && quarto render publications.qmd && quarto render index.qmd && cd ..

# 4. Review, commit and push the site branch.
cd docs/_site
git status --short     # expect: the changed pages, search.json, sitemap.xml (lastmod only)
git add -A && git commit -m "Deploy docs: <what changed>" && git push origin gh-pages
```

`docs/_site/` is gitignored on `master`, so the nested `gh-pages` clone never leaks into source commits.
Do not hand-edit generated HTML; fix the source and re-render. If you instead render on a machine that
already has a full local `docs/_site` (with `_freeze/`), you can copy it into a `gh-pages` clone with
`rsync -a --delete --exclude .git --exclude .DS_Store docs/_site/ <gh-pages-clone>/`, then commit and push.

### After pushing

1. Pages builds automatically (workflow "pages build and deployment"): `gh run list -R mofography/mofography -L 1`.
   It usually takes a few minutes.
2. If the deploy job fails with `Deployment failed, try again later`, rerun it: `gh run rerun <id>`.
   If that fails too, push an empty commit to `gh-pages`
   (`git commit --allow-empty -m "Retrigger Pages deploy" && git push`). The previous site keeps being
   served meanwhile.
3. Verify on the live URL, e.g.
   `curl -s https://mofography.github.io/mofography/publications.html | grep "<new text>"`.

### Rules and gotchas

- `gh-pages` needs a `.nojekyll` file at its root (already present). Keep it.
- `docs/_quarto.yml` renders `*.qmd`, `*.ipynb` and `*.md` in `docs/`. Files that must not be published are
  excluded there with `"!file.md"` (e.g. `PAPER_METADATA_TODO.md`, the maintainer checklist). Add new
  internal notes to that list or keep them outside `docs/`.
- `site-url` in `docs/_quarto.yml` is `https://mofography.github.io/mofography/` (org `mofography`, repo `mofography`).
- The PoroNet tutorial downloads `Pore_Graph_for_H2.pkl` from the release `poronet-data-v1`. Keep that
  release and tag intact; do not commit the file.
- Do not commit rendered output (`_site`, `_build`, `_freeze`, `.jupyter_cache`, `outputs/`) or large generated
  HTML/binaries to `master`. The history was rewritten once (Nov 2026) to remove them.
- Citation changes go in four places together: `README.md`, `CITATION.cff`, `docs/index.qmd` (Cite section)
  and `docs/publications.qmd`. Methods-paper author order: Arun Gopalan, Chao Zheng, Behnam Azadpour,
  Kaihang Shi, Lev Sarkisov.
