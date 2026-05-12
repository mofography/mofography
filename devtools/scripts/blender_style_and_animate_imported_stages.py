"""
Style and animate manually imported MOFOGRAPHY stage GLBs in Blender.

How to use (UI):
1) Open your Blender file where you already imported all stage GLBs.
2) Open Scripting tab.
3) New text block -> paste this script (or open this file).
4) Click "Run Script".

The script expects object names that include stage prefixes like:
  01_framework, 02_cell, 03_regions, ..., 09_maxima_pbc
which is how the GLB exporter names objects.
"""

from __future__ import annotations

import bpy
from mathutils import Vector


# -------------------------
# User-tweakable settings
# -------------------------
FPS = 30
HOLD_FRAMES_PER_STAGE = 28
FADE_FRAMES = 12
FINAL_HOLD = 30
ROTATIONS = 1.0  # full turns around Z
CAMERA_LOCATION = (35.0, -35.0, 22.0)
SUN_LOCATION = (10.0, -8.0, 16.0)
SUN_ENERGY = 4.0

STAGE_ORDER = [
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

# Simple, readable colors for explanation animation
STAGE_BASE_COLORS = {
    "01_framework": (0.78, 0.78, 0.78, 1.0),      # light gray
    "02_cell": (1.00, 0.45, 0.08, 1.0),           # orange
    "03_regions": (0.45, 0.70, 1.00, 1.0),        # blue-ish
    "04_regions_pbc": (0.86, 0.55, 1.00, 1.0),    # violet
    "05_pocket_types": (0.40, 0.90, 0.55, 1.0),   # green
    "06_connections": (1.00, 0.20, 0.20, 1.0),    # red
    "07_maxima_distinct": (1.00, 0.95, 0.15, 1.0),# yellow
    "08_windows": (0.20, 1.00, 1.00, 1.0),        # cyan
    "09_maxima_pbc": (1.00, 0.35, 0.80, 1.0),     # magenta
}

# Slightly transparent surfaces for readability
STAGE_ALPHA = {
    "03_regions": 0.35,
    "04_regions_pbc": 0.25,
    "05_pocket_types": 0.35,
}

STAGES_USE_VERTEX_COLORS = {"01_framework", "03_regions", "04_regions_pbc", "05_pocket_types"}


def get_mesh_objects():
    return [o for o in bpy.data.objects if o.type == "MESH"]


def stage_for_object(obj_name: str) -> str | None:
    for stage in STAGE_ORDER:
        if stage in obj_name:
            return stage
    return None


def ensure_material(obj, stage: str):
    if obj.data.materials and obj.data.materials[0] is not None:
        mat = obj.data.materials[0]
    else:
        mat = bpy.data.materials.new(name=f"{stage}_mat")
        obj.data.materials.append(mat)

    mat.use_nodes = True
    mat.blend_method = "BLEND"
    if hasattr(mat, "shadow_method"):
        mat.shadow_method = "HASHED"

    nt = mat.node_tree
    nodes = nt.nodes
    links = nt.links
    bsdf = nodes.get("Principled BSDF")
    if bsdf is None:
        bsdf = nodes.new("ShaderNodeBsdfPrincipled")
    output = nodes.get("Material Output")
    if output is None:
        output = nodes.new("ShaderNodeOutputMaterial")
    if not bsdf.outputs["BSDF"].is_linked:
        links.new(bsdf.outputs["BSDF"], output.inputs["Surface"])

    base = STAGE_BASE_COLORS.get(stage, (0.8, 0.8, 0.8, 1.0))
    alpha = STAGE_ALPHA.get(stage, 1.0)

    # Preserve vertex colors from exported meshes where meaningful:
    # - framework (CPK atoms + bond colors from mesh_for_atoms)
    # - three region-derived isosurfaces (region labels -> colorscale mapping)
    has_color_attr = len(obj.data.color_attributes) > 0
    if stage in STAGES_USE_VERTEX_COLORS and has_color_attr:
        vcol = nodes.get("VertexColor")
        if vcol is None:
            vcol = nodes.new("ShaderNodeVertexColor")
            vcol.name = "VertexColor"
        # Prefer attribute named "Color" (seen in GLB import), fallback to first one.
        attr_name = "Color"
        if attr_name not in obj.data.color_attributes:
            attr_name = obj.data.color_attributes[0].name
        vcol.layer_name = attr_name
        links.new(vcol.outputs["Color"], bsdf.inputs["Base Color"])
        if stage == "01_framework":
            # Framework stays opaque; CPK + bond coloring is in vertex colors.
            alpha = 1.0
    else:
        bsdf.inputs["Base Color"].default_value = base

    bsdf.inputs["Roughness"].default_value = 0.35
    bsdf.inputs["Specular IOR Level"].default_value = 0.2
    bsdf.inputs["Alpha"].default_value = alpha
    return mat, bsdf, alpha


def clear_existing_animation():
    for obj in bpy.data.objects:
        if obj.animation_data:
            obj.animation_data_clear()
    for mat in bpy.data.materials:
        if mat.animation_data:
            mat.animation_data_clear()


def create_camera_and_pivot(target: Vector):
    cam = bpy.data.objects.get("Camera_Main")
    if cam is None:
        bpy.ops.object.camera_add(location=CAMERA_LOCATION)
        cam = bpy.context.active_object
        cam.name = "Camera_Main"

    pivot = bpy.data.objects.get("Camera_Pivot")
    if pivot is None:
        bpy.ops.object.empty_add(type="PLAIN_AXES", location=target)
        pivot = bpy.context.active_object
        pivot.name = "Camera_Pivot"
    else:
        pivot.location = target

    cam.parent = pivot
    cam.matrix_parent_inverse = pivot.matrix_world.inverted()

    track = None
    for c in cam.constraints:
        if c.type == "TRACK_TO":
            track = c
            break
    if track is None:
        track = cam.constraints.new(type="TRACK_TO")
    track.target = pivot
    track.track_axis = "TRACK_NEGATIVE_Z"
    track.up_axis = "UP_Y"
    return cam, pivot


def ensure_sun():
    sun = bpy.data.objects.get("Sun_Key")
    if sun is None:
        bpy.ops.object.light_add(type="SUN", location=SUN_LOCATION)
        sun = bpy.context.active_object
        sun.name = "Sun_Key"
    sun.location = SUN_LOCATION
    sun.data.energy = SUN_ENERGY
    return sun


def keyframe_stage_reveal(stage_objects: dict[str, list[bpy.types.Object]]):
    frame = 1
    for stage in STAGE_ORDER:
        objs = stage_objects.get(stage, [])
        for obj in objs:
            obj.hide_viewport = True
            obj.hide_render = True
            obj.keyframe_insert(data_path="hide_viewport", frame=1)
            obj.keyframe_insert(data_path="hide_render", frame=1)

        for obj in objs:
            obj.hide_viewport = False
            obj.hide_render = False
            obj.keyframe_insert(data_path="hide_viewport", frame=frame)
            obj.keyframe_insert(data_path="hide_render", frame=frame)

            stage = stage_for_object(obj.name) or "01_framework"
            _, bsdf, target_alpha = ensure_material(obj, stage)
            bsdf.inputs["Alpha"].default_value = 0.0
            bsdf.inputs["Alpha"].keyframe_insert(data_path="default_value", frame=frame)
            bsdf.inputs["Alpha"].default_value = target_alpha
            bsdf.inputs["Alpha"].keyframe_insert(data_path="default_value", frame=frame + FADE_FRAMES)

        frame += HOLD_FRAMES_PER_STAGE

    return frame + FINAL_HOLD


def linearize_fcurve_keys(obj):
    if obj.animation_data and obj.animation_data.action:
        for fc in obj.animation_data.action.fcurves:
            for kp in fc.keyframe_points:
                kp.interpolation = "LINEAR"


def main():
    mesh_objects = get_mesh_objects()
    if not mesh_objects:
        raise RuntimeError("No mesh objects found. Import GLBs first, then run this script.")

    clear_existing_animation()

    # Group objects by stage prefix found in names
    stage_objects: dict[str, list[bpy.types.Object]] = {k: [] for k in STAGE_ORDER}
    for obj in mesh_objects:
        stage = stage_for_object(obj.name)
        if stage is not None:
            stage_objects[stage].append(obj)

    # Warn if any stage is empty
    for stage in STAGE_ORDER:
        if not stage_objects[stage]:
            print(f"[WARN] No objects found for stage: {stage}")

    # Prepare base materials first
    for stage, objs in stage_objects.items():
        for obj in objs:
            ensure_material(obj, stage)

    # Camera target from object centroid
    center = sum((o.location for o in mesh_objects), Vector((0.0, 0.0, 0.0))) / max(len(mesh_objects), 1)
    cam, pivot = create_camera_and_pivot(center)
    ensure_sun()

    end_frame = keyframe_stage_reveal(stage_objects)

    # Camera orbit
    pivot.rotation_euler[2] = 0.0
    pivot.keyframe_insert(data_path="rotation_euler", frame=1)
    pivot.rotation_euler[2] = 6.283185307 * ROTATIONS
    pivot.keyframe_insert(data_path="rotation_euler", frame=end_frame)
    linearize_fcurve_keys(pivot)

    # Scene render settings
    scene = bpy.context.scene
    scene.camera = cam
    scene.render.fps = FPS
    scene.frame_start = 1
    scene.frame_end = end_frame
    scene.render.resolution_x = 1920
    scene.render.resolution_y = 1080

    engines = {item.identifier for item in scene.render.bl_rna.properties["engine"].enum_items}
    if "BLENDER_EEVEE_NEXT" in engines:
        scene.render.engine = "BLENDER_EEVEE_NEXT"
    elif "BLENDER_EEVEE" in engines:
        scene.render.engine = "BLENDER_EEVEE"
    else:
        scene.render.engine = "CYCLES"

    scene.render.image_settings.file_format = "FFMPEG"
    scene.render.ffmpeg.format = "MPEG4"
    scene.render.ffmpeg.codec = "H264"
    scene.render.ffmpeg.constant_rate_factor = "HIGH"
    scene.render.ffmpeg.ffmpeg_preset = "GOOD"

    print("Done. Styled objects, keyframed stage reveal, and configured render settings.")


if __name__ == "__main__":
    main()

