"""Group the quartodoc reference by status (docs/function_status.yml) inside each module group.

Pipeline (repo root):
    python scripts/dev/merge_function_status.py        # proposals + overrides -> docs/function_status.yml
    python scripts/dev/split_deprecated_reference.py   # rewrites the quartodoc block of docs/_quarto.yml
    cd docs && quartodoc build && python ../scripts/dev/tag_reference.py

Idempotent: earlier sub-sections are merged back into their module first. Only the ``quartodoc:``
block at the end of docs/_quarto.yml is rewritten (it must be the last top-level key).
"""
from pathlib import Path

import yaml

DOCS = Path(__file__).resolve().parents[2] / "docs"
QUARTO = DOCS / "_quarto.yml"
ORDER = [  # (status, sub-section title, description)
    ("recommended", "Recommended", "Go-to functions for this module."),
    ("alternate", "Alternate", "Valid, but not the default choice."),
    ("utility", "Other", "Low-level helpers and building blocks."),
    ("experimental", "Experimental", "Prototypes; may change or be removed."),
    ("legacy", "Legacy", "Superseded by a better function; still works."),
    ("deprecated", "Deprecated", "Do not use in new code."),
]
RENAME = {"Visualize": "Plotly visualization", "Three.js viz": "Three.js visualization"}


def main():
    text = QUARTO.read_text()
    head, _, _ = text.partition("\nquartodoc:\n")
    cfg = yaml.safe_load(text)["quartodoc"]
    status = yaml.safe_load((DOCS / "function_status.yml").read_text())

    modules = []
    for s in cfg["sections"]:
        if "title" in s:
            s.pop("contents_all", None)
            modules.append(s)
        else:
            modules[-1]["contents"] += s["contents"]

    out = []
    for s in modules:
        pkg, names = s["package"], [n if isinstance(n, str) else n["name"] for n in s["contents"]]
        s["title"] = RENAME.get(s["title"], s["title"])
        s["contents"] = []
        out.append(s)
        for key, title, desc in ORDER:
            group = [n for n in names if status[f"{pkg}.{n}"]["status"] == key]
            if group:
                out.append({"subtitle": title, "desc": desc, "package": pkg, "contents": group})
    cfg["sections"] = out
    QUARTO.write_text(head + "\n" + yaml.safe_dump({"quartodoc": cfg}, sort_keys=False, allow_unicode=True, width=100))
    print(f"{len(out)} sections written")


if __name__ == "__main__":
    main()
