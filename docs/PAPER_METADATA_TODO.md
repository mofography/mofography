# Paper & metadata TODO (maintainers)

Complete these when the manuscript is on arXiv / accepted and when the org repo goes public.

- [ ] **arXiv:** Replace `PLACEHOLDER_ARXIV_ID` in `README.md` and `CITATION.cff` with the real ID and URL.
- [ ] **Journal:** Add volume, pages, and publisher DOI when accepted (update `preferred-citation` in `CITATION.cff`).
- [ ] **Zenodo:** After `mofography/mofography` is public, follow [releases-and-citing.md](releases-and-citing.md) → Zenodo checklist; add version DOI to `CITATION.cff` and README.
- [ ] **CITATION.cff authors:** Software entry is currently minimal; align `authors` / `preferred-citation.authors` with the final author list from the paper.
- [ ] **Badges:** Refresh CI / codecov URLs if workflows move fully to the org repo.
- [ ] **Install:** If the project adopts a root `pyproject.toml` or publishes wheels, update README install paths and remove duplicate instructions in `docs/install.qmd`.
- [ ] **mgr_cli_py:** When [arung-northwestern/mgr_cli_py](https://github.com/arung-northwestern/mgr_cli_py) is public, confirm README links and add PyPI/docker URLs if published.
- [ ] **HF demo:** Confirm Space URL and resource limits in README if hosting changes.
- [ ] **GitHub Pages:** Deploy Quarto output (`cd docs && quarto render`) to Pages for `mofography.github.io/mofography/`; enable badge link in README.

Source manuscript (private): `mgr_sn` repo, `mgr_sn_v2.org`.
