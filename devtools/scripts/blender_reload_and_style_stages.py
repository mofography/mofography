"""
Reload staged GLBs and apply styling/animation in the current Blender file.

Run this inside Blender's Scripting workspace.
"""

from __future__ import annotations

from pathlib import Path
import runpy
import bpy


STAGE_DIR = Path("/Users/arungopalan/repos/mofography/outputs/blender_stages")
STAGE_PREFIXES = [
    "01_framework",
    "02_cell",
    "03_regions",
    "04_regions_pbc",
    "05_pocket_types",
    "06_connections",
    "07_maxima_distinct",
    "08_windows",
    "09_maxima_pbc",
]


def main():
    # Remove existing stage objects only
    for obj in list(bpy.data.objects):
        if any(prefix in obj.name for prefix in STAGE_PREFIXES):
            bpy.data.objects.remove(obj, do_unlink=True)

    # Import newest GLBs
    for glb in sorted(STAGE_DIR.glob("*.glb")):
        bpy.ops.import_scene.gltf(filepath=str(glb))
        print(f"Imported {glb.name}")

    # Apply materials/colors/transparency/animation
    runpy.run_path(
        "/Users/arungopalan/repos/mofography/devtools/scripts/blender_style_and_animate_imported_stages.py",
        run_name="__main__",
    )
    print("Reload + style complete.")


if __name__ == "__main__":
    main()

