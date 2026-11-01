"""Post-process `quartodoc build` output: status badges, filter bar and top note.

Run from docs/ after `quartodoc build`:  python ../scripts/dev/tag_reference.py
Reads docs/function_status.yml (see scripts/dev/merge_function_status.py).
"""
import re
from pathlib import Path

import yaml

REF = Path("reference")
STATUS = yaml.safe_load(Path("function_status.yml").read_text())
LABEL = {"recommended": "go-to", "alternate": "alternate", "utility": "other",
         "experimental": "experimental", "legacy": "legacy", "deprecated": "deprecated"}
ORDER = ["recommended", "alternate", "utility", "experimental", "legacy", "deprecated"]

NOTE = """::: {.callout-note appearance="simple"}
You probably won't need most of these functions. The core workflow is covered in the
[Learn tutorials](../learn.qmd); this reference lists everything the library exposes.
Some functions are **experimental**, **legacy** or **deprecated** and are marked as such.
Use the tags below to filter, or the groups inside each module.
:::

"""


def badge(status):
    return f'<span class="mg-tag mg-tag-{status}">{LABEL[status]}</span>'


def filter_bar():
    btns = ['<button class="mg-tag-btn active" data-tag="all">all</button>']
    btns += [f'<button class="mg-tag-btn mg-tag-{s}" data-tag="{s}">{LABEL[s]}</button>' for s in ORDER]
    return '```{=html}\n<div class="mg-tag-filter">Show: ' + " ".join(btns) + "</div>\n```\n\n"


def tag_index():
    p = REF / "index.qmd"
    text = p.read_text()
    if "mg-tag-filter" in text:
        return
    row = re.compile(r"^\| (\[[^\]]+\]\([^#)]+#([\w.]+)\)) \|", re.M)

    def sub(m):
        key = m.group(2).rsplit(".", 1)
        s = STATUS.get(".".join(key))
        return f"| {m.group(1)} {badge(s['status'])} |" if s else m.group(0)

    text = row.sub(sub, text)
    text = re.sub(r":(?:func|class|meth|mod|attr):`~?([^`]+)`", r"`\1`", text)  # strip Sphinx roles
    head, sep, rest = text.partition("\n")
    p.write_text(head + "\n\n" + NOTE + filter_bar() + rest)


def tag_pages():
    for p in REF.glob("*.qmd"):
        if p.name == "index.qmd":
            continue
        text = p.read_text()
        m = re.match(r"# .+? \{ #([\w.]+) \}\n", text)
        if not m or "mg-tag" in text:
            continue
        s = STATUS.get(m.group(1).rsplit(".", 1)[0] + "." + m.group(1).rsplit(".", 1)[1])
        if not s or s["status"] == "utility":
            continue
        line = badge(s["status"])
        if s.get("note"):
            line += f" {s['note']}"
        if s.get("use"):
            line += f" Use `{s['use']}` instead." if s["status"] in ("legacy", "deprecated", "alternate") else ""
        p.write_text(text[: m.end()] + f"\n{line}\n" + text[m.end():])


if __name__ == "__main__":
    tag_index()
    tag_pages()
    print("tagged")
