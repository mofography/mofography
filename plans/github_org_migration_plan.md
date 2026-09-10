# GitHub org migration plan — personal sandbox + private org copy

> **Strategy:** Keep `arung-northwestern/mofography` (private sandbox). Seed new
> `mofography/mofography` (private org copy) from a **curated `master`**. Full
> commit history preserved — not squashed. Docs website polish **after** org seed.
>
> **Status:** Groundwork / planning (2026-09-10). Do not create org until Phase 1
> cleanup is committed and pushed to personal remote.

---

## Goals

| Goal | How |
|------|-----|
| Personal repo stays | No transfer. `origin` = personal forever (or until you choose otherwise) |
| Org repo for group / future public | Second remote `org`; seed once, promote with `git push org master` |
| Full commit history + lead-dev attribution | Normal `git push` — every commit, author, date preserved |
| Exclude dev clutter from org | Untrack bloat on `master` **before** first org push |
| arXiv gate | Both stay private now; flip org public later |
| Existing collaborator forks | Stay on personal repo until you migrate them (post-arXiv) |

---

## Current repo audit (2026-09-10)

### Git state

| Item | Value |
|------|-------|
| Remote | `origin` → `arung-northwestern/mofography` (private) |
| Branch | `master` (17 commits **unpushed** to origin) |
| Tip | `e60e67b` |
| Total commits | ~643 |
| Tags | none |
| GitHub fork network | 3 private forks: `Andrew-S-Rosen`, `Shi-Research-Group`, `ConorNCleeton` |
| `.git` size | ~3 GB (historical blobs: blend files, large HTML movies, Playwright browsers) |
| Commit author | `arungopalan2020@u.northwestern.edu` (verify linked on GitHub profile) |

### Tracked file counts by top-level dir

| Dir | Tracked files | Verdict for org seed |
|-----|--------------:|----------------------|
| `outputs/` | **3542** (~926 MB, mostly PNG animation frames) | **DROP from git** |
| `docs/` | 646 (~1 GB with caches) | **KEEP source**, drop caches |
| `examples/` | 79 | **KEEP** scripts + CIFs; drop heavy viz HTML |
| `internal_docs/` | 61 | **DROP from git** (dev notes, Shi handoff) |
| `mofography/` | 56 | **KEEP** |
| `scripts/` | 30 | **KEEP** |
| `devtools/` | 23 | **DROP from git** (animation tooling) |
| `mgr_rust/` | 13 | **KEEP** |
| `plans/` | 9 | **DROP from git** (bug investigations) |
| `tests/` | 6 | **KEEP** |
| `iframe_figures/` | 2 (~35 MB HTML) | **DROP from git** |
| `.cursor/` | 2 | **DROP from git** (agent config — personal) |
| `parity_outputs/` | 1 | **DROP from git** |

### Largest tracked files (working tree)

| File | Size | Action |
|------|-----:|--------|
| `docs/framework_visualization_example.ipynb` | 64 MB | Keep for now (under 100 MB limit); strip outputs later per website plan |
| `docs/.jupyter_cache/**` | 42 MB × many | Untrack |
| `docs/_freeze/**` | 21 MB × many | Untrack (Quarto regen on render) |
| `mofography/examples/irmof_example*.ipynb` | 31–34 MB | Consider untrack duplicates (also in `docs/`) |
| `iframe_figures/*.html` | 17 MB each | Untrack |
| `outputs/**` | 926 MB total | Untrack entire tree |

### Already untracked (good — stays local only)

- `examples/viz_outputs/` (most HTML/PDF/NPZ outputs)
- `parity_outputs/` (except 1 tracked file)
- `plans/` (newer plans: website rebuild, tutorial migration, etc.)
- `.tmp_*`, `.doom-cockpit/`, `.omp/`
- `graphify-out/`, `.code-review-graph/` (gitignored)

### Docs readiness (defer to post-org)

- Quarto config exists (`docs/_quarto.yml`, freeze mode on)
- Site not publication-ready (website rebuild plan Phase 1–4 open)
- **OK to seed org with docs source as-is** — private org, no public site yet
- Stale `readthedocs.yml` (Sphinx-era) → remove or ignore
- README badges still say `REPLACE_WITH_OWNER_ACCOUNT`

---

## What goes in org repo vs personal-only

### Org seed includes (core product)

```
mofography/          # Python library
mgr_rust/            # Rust extensions
tests/
scripts/             # incl. scripts/dev/ — useful for repro
examples/            # scripts + TestCIFs (not viz_outputs HTML dumps)
docs/                # Quarto source (notebooks, qmd, reference) — no _site/cache/freeze
.github/             # CONTRIBUTING, CODE_OF_CONDUCT, future CI
LICENSE, README.md, setup.py, setup.cfg, requirements.txt, pytest.ini, MANIFEST.in
```

### Exclude from org (and ideally untrack on master)

| Path | Why |
|------|-----|
| `outputs/` | 900+ PNG frames per animation; reproducible from scripts |
| `docs/.jupyter_cache/`, `docs/_freeze/` | Build artifacts; Quarto regenerates |
| `docs/_site/` | Already gitignored |
| `internal_docs/` | **Not junk** — valuable dev/research notes, handoffs, one-off test scripts, Shi fork changelog. Wrong audience for org product repo; keep locally. Some content will be promoted to `docs/` later (see website rebuild plan). |
| `plans/` | Bug reports, watershed debug, agent plans |
| `devtools/` | Blender/animation dev scripts |
| `iframe_figures/` | Large standalone HTML embeds |
| `parity_outputs/` | Dev parity HTML |
| `.cursor/`, `AGENTS.md` | Agent/tooling config — no secrets in tracked files (verified), but personal workflow |
| `readthedocs.yml` | Stale Sphinx config |
| `core2000_window_fp_accessible_comparison.csv` | One-off benchmark artifact |

### Privacy & secrets (verified 2026-09-10)

- Org repo created with `--private`; stays private until you manually change visibility post-arXiv
- `git push org master` does **not** change visibility
- Scanned tracked files: no API keys, tokens, passphrases, `.env` files
- Already gitignored: `.letta/`, `.serena/`, `.claude/`, `.cursor/mcp.json`, `graphify-out/`
- Untracked local only: `.doom-cockpit/`, `.omp/`, `.tmp_*`
- After seed: verify org repo Settings → Danger Zone shows **Private**

### Personal-only after cleanup

Anything gitignored stays on disk locally but never reaches org:

- Local viz HTML under `examples/viz_outputs/`
- `parity_outputs/`, `.tmp_*`, agent tooling
- Extra plans in untracked `plans/`
- `graphify-out/` (regenerate locally)

**Important:** Personal and org share the same `master` if you use dual-remote push. "Personal-only" = **gitignored local files**, not separate tracked trees. For tracked experiments you want on personal but not org, use a `sandbox` branch on personal only (advanced — not needed if junk is untracked).

---

## Phase checklist

### Phase 0 — You verify in browser (5 min)

- [ ] GitHub profile has `arungopalan2020@u.northwestern.edu` in emails (commit attribution)
- [ ] Confirm `arung-northwestern/mofography` visibility = Private
- [ ] Pick org name: `mofography` (appears available)
- [ ] Decide org members now or later (can be solo initially)

### Phase 1 — Repo cleanup on personal (local + commit)

**Goal:** Slim `master` before any org seed. One focused commit.

1. **Extend `.gitignore`** — add:
   ```
   outputs/
   docs/.jupyter_cache/
   docs/_freeze/
   internal_docs/
   plans/
   devtools/
   iframe_figures/
   parity_outputs/
   .cursor/
   examples/viz_outputs/
   archived_tutorials/
   .tmp_*
   .doom-cockpit/
   .omp/
   ```

2. **Untrack (keep files on disk):**
   ```bash
   git rm -r --cached outputs/ docs/.jupyter_cache docs/_freeze \
     internal_docs/ plans/ devtools/ iframe_figures/ parity_outputs/ .cursor/
   ```

3. **Optional same commit:**
   - Remove stale `readthedocs.yml`
   - Fix README clone URL (still personal for now; org URL after seed)
   - Stage unstaged library work (`mofography/rust_fun.py`, etc.) if intended

4. **Commit:** `chore: untrack dev artifacts before org seed`

5. **Push personal:** `git push origin master` (includes 17 pending feature commits)

**Does NOT rewrite history.** Old blobs still in `.git` (~3 GB). Optional Phase 1b later: `git filter-repo` to shrink clones before org seed.

### Phase 2 — Create org + empty repo (browser or gh cli)

```bash
# Browser: github.com/organizations/plan → create free org "mofography"

# Or CLI (after org exists):
gh org create mofography --description "MOF pore analysis and visualization"  # if supported
gh repo create mofography/mofography --private --description "MOFography library"
```

You do org billing/plan confirmation in browser once.

### Phase 3 — Seed org (one-time)

```bash
git remote add org https://github.com/mofography/mofography.git
git push org master
# optional: git push org --all   # only if you need other branches on org
```

Verify in browser:
- [ ] Commit count ~643+ (not 1)
- [ ] Your author on recent commits
- [ ] No `outputs/` tree in org repo file browser
- [ ] Repo stays private

### Phase 4 — Dual-remote workflow (ongoing)

```bash
git push origin master    # always — personal sandbox
git push org master       # when ready — promote to org
```

Optional git config:
```bash
git config remote.pushDefault origin   # default push = personal only
```

### Phase 5 — Post-seed (can wait)

- [ ] Update README / `docs/_quarto.yml` GitHub links → `mofography/mofography`
- [ ] Rewrite `.github/CONTRIBUTING.md` for org (not generic fork boilerplate)
- [ ] Website rebuild plan (Quarto site, CI, Pages) — when ready
- [ ] arXiv drop → flip **org** repo public (personal stays private)
- [ ] Tell Shi / Rosen / Cleeton forks to re-target org repo when public

---

## History bloat note

Untracking `outputs/` removes **future** weight. Past commits still contain:
- 290 MB `.blend` file
- Playwright browser binaries
- 135 MB HTML movie files

Org seed copies that history too. Acceptable for private org short-term. To fix before seed:

```bash
# Phase 1b (optional, destructive to SHAs on rewritten commits)
# pip install git-filter-repo
# git filter-repo --path outputs/ --invert-paths
# git push origin master --force   # coordinate if others clone personal repo
```

Only do 1b if clone size matters before inviting collaborators to org.

---

## What we are NOT doing

| Action | Why not |
|--------|---------|
| Transfer personal → org | Loses personal sandbox |
| GitHub "Import repository" | Fresh repo, awkward history |
| Squash / orphan branch | Destroys lead-dev commit trail |
| Public org repo now | arXiv gate |
| Docs website before seed | You said OK to defer |
| Touch existing collaborator forks | They stay on personal until you migrate |

---

## Execution order (recommended)

```
Phase 0 (browser verify)
    ↓
Phase 1 (cleanup commit + push personal)  ← do this first
    ↓
Phase 2 (create org + empty repo)
    ↓
Phase 3 (seed org, verify history)
    ↓
Phase 4 (dual-remote habit)
    ↓
Phase 5 (docs, public flip) — later
```

---

## Quick reference — sizes removed from org tip if Phase 1 done

| Removed from tip | ~Size |
|------------------|------:|
| `outputs/` | 926 MB |
| `docs/.jupyter_cache` | 268 MB |
| `docs/_freeze` | 154 MB |
| `internal_docs/` | 21 MB |
| `iframe_figures/` | 35 MB |
| **Total off tip** | **~1.4 GB** |

Tracked files drop from ~4488 → ~900 (mostly `docs/` + library).
