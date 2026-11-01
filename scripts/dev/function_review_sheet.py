"""Write a review sheet of every documented function so you can mark keep / legacy / remove.

    python scripts/dev/function_review_sheet.py > internal_docs/docs_redesign/function_review.md

Put an x in one column per row, then hand the file back (or tell the agent which rows to mark).
Marking a function "Legacy" means: first docstring line becomes ``Legacy: use <x>.`` and
``split_deprecated_reference.py`` moves it to the bottom block of its module group.
"""
import importlib
import inspect
import warnings
from pathlib import Path

import yaml

warnings.filterwarnings("ignore")
cfg = yaml.safe_load((Path(__file__).resolve().parents[2] / "docs" / "_quarto.yml").read_text())["quartodoc"]
print("# Function review sheet\n\nMark one column per row: **keep**, **legacy** (old workflow, confusing), **remove**.\n")
seen = {}
for s in cfg["sections"]:
    if "title" in s:
        title = s["title"]
        seen[title] = s
        print(f"\n## {title}\n\n| function | first docstring line | keep | legacy | remove |\n|---|---|:-:|:-:|:-:|")
    mod = importlib.import_module(s["package"])
    for n in s["contents"]:
        n = n if isinstance(n, str) else n["name"]
        doc = (inspect.getdoc(getattr(mod, n, None)) or "").strip().split("\n")[0][:90].replace("|", "/")
        dep = " (deprecated)" if "subtitle" in s else ""
        print(f"| `{n}`{dep} | {doc} | [ ] | [ ] | [ ] |")
