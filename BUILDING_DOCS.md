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

## Publishing on GitHub Pages (once the repo is public)

1. Set `site-url` in `docs/_quarto.yml` to the final address, for example `https://<org>.github.io/<repo>/`.
2. Render, then publish the rendered site to the `gh-pages` branch:

   ```bash
   cd docs
   quarto publish gh-pages
   ```

3. In the repository settings, under **Pages**, choose the `gh-pages` branch as the source.
4. To update the site later, run `quarto render` and `quarto publish gh-pages` again.

This publishes from your machine, so no CI workflow is needed.
