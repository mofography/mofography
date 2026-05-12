"""
Run inside Blender to import staged GLBs and make a reveal animation.

CLI example:
  blender --python devtools/scripts/blender_import_staged_glb_and_animate.py -- \
      --stages outputs/blender_stages --output outputs/blender_stages/tutorial_reveal.blend
"""

from __future__ import annotations

import argparse
from pathlib import Path

import bpy
from mathutils import Vector


def _parse_args():
    import sys

    argv = sys.argv
    if "--" in argv:
        argv = argv[argv.index("--") + 1 :]
    else:
        argv = []
    parser = argparse.ArgumentParser()
    parser.add_argument("--stages", type=Path, default=Path("outputs/blender_stages"))
    parser.add_argument("--output", type=Path, default=Path("outputs/blender_stages/tutorial_reveal.blend"))
    parser.add_argument("--fps", type=int, default=30)
    parser.add_argument("--hold", type=int, default=25, help="Frames to hold each stage before next reveal")
    parser.add_argument("--fade", type=int, default=10, help="Frames for transparency transition")
    return parser.parse_args(argv)


def _clear_scene():
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete(use_global=False)
    for block in bpy.data.meshes:
        if block.users == 0:
            bpy.data.meshes.remove(block)
    for block in bpy.data.materials:
        if block.users == 0:
            bpy.data.materials.remove(block)


def _import_stage(path: Path):
    before = set(bpy.data.objects.keys())
    bpy.ops.import_scene.gltf(filepath=str(path))
    after = set(bpy.data.objects.keys())
    new_names = after - before
    return [bpy.data.objects[name] for name in new_names if bpy.data.objects[name].type == "MESH"]


def _ensure_alpha_material(obj):
    if not obj.data.materials:
        mat = bpy.data.materials.new(name=f"{obj.name}_mat")
        mat.use_nodes = True
        obj.data.materials.append(mat)
    mat = obj.data.materials[0]
    mat.use_nodes = True
    mat.blend_method = "BLEND"
    if hasattr(mat, "shadow_method"):
        mat.shadow_method = "HASHED"

    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf is None:
        return None
    return bsdf


def _keyframe_opacity(obj, frame_a: int, alpha_a: float, frame_b: int, alpha_b: float):
    bsdf = _ensure_alpha_material(obj)
    if bsdf is None:
        return
    bsdf.inputs["Alpha"].default_value = alpha_a
    bsdf.inputs["Alpha"].keyframe_insert(data_path="default_value", frame=frame_a)
    bsdf.inputs["Alpha"].default_value = alpha_b
    bsdf.inputs["Alpha"].keyframe_insert(data_path="default_value", frame=frame_b)


def _create_camera_rig(target: Vector):
    bpy.ops.object.camera_add(location=(35, -35, 22))
    cam = bpy.context.active_object
    cam.name = "Camera_Main"

    bpy.ops.object.empty_add(type="PLAIN_AXES", location=target)
    pivot = bpy.context.active_object
    pivot.name = "Camera_Pivot"
    cam.parent = pivot
    cam.matrix_parent_inverse = pivot.matrix_world.inverted()

    constraint = cam.constraints.new(type="TRACK_TO")
    constraint.target = pivot
    constraint.track_axis = "TRACK_NEGATIVE_Z"
    constraint.up_axis = "UP_Y"
    return cam, pivot


def _set_render_defaults(fps: int, end_frame: int):
    scene = bpy.context.scene
    available = {item.identifier for item in scene.render.bl_rna.properties["engine"].enum_items}
    if "BLENDER_EEVEE_NEXT" in available:
        scene.render.engine = "BLENDER_EEVEE_NEXT"
    elif "BLENDER_EEVEE" in available:
        scene.render.engine = "BLENDER_EEVEE"
    else:
        scene.render.engine = "CYCLES"
    scene.render.fps = fps
    scene.frame_start = 1
    scene.frame_end = end_frame


def main():
    args = _parse_args()
    stage_dir = args.stages.resolve()
    stage_files = sorted(stage_dir.glob("*.glb"))
    if not stage_files:
        raise RuntimeError(f"No .glb files found in {stage_dir}")

    _clear_scene()

    stage_objects = {}
    for f in stage_files:
        objs = _import_stage(f)
        stage_objects[f.stem] = objs
        for o in objs:
            o.hide_render = True
            o.hide_viewport = True
            o.keyframe_insert(data_path="hide_render", frame=1)
            o.keyframe_insert(data_path="hide_viewport", frame=1)

    all_meshes = [o for objs in stage_objects.values() for o in objs]
    center = Vector((0.0, 0.0, 0.0))
    if all_meshes:
        center = sum((o.location for o in all_meshes), Vector((0, 0, 0))) / len(all_meshes)

    _, pivot = _create_camera_rig(center)

    frame = 1
    for _, objs in stage_objects.items():
        for o in objs:
            o.hide_render = False
            o.hide_viewport = False
            o.keyframe_insert(data_path="hide_render", frame=frame)
            o.keyframe_insert(data_path="hide_viewport", frame=frame)
            _keyframe_opacity(o, frame, 0.0, frame + args.fade, 1.0)
        frame += args.hold

    pivot.rotation_euler[2] = 0.0
    pivot.keyframe_insert(data_path="rotation_euler", frame=1)
    pivot.rotation_euler[2] = 6.283185307
    pivot.keyframe_insert(data_path="rotation_euler", frame=frame)

    _set_render_defaults(args.fps, frame + 15)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.save_as_mainfile(filepath=str(args.output))
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
