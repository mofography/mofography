#!/usr/bin/env python3
"""
Create an interactive Three.js HTML viewer for framework + segmentation + MD snapshots.

Current composition:
- Framework atoms (from CIF) as one trace
- Segmentation regions surface (from CIF dgrid/regions) as one trace
- Guest atoms per MD snapshot (from RASPA-style PDB movie) as one animated trace

Controls:
- Trace visibility toggles
- Snapshot slider + play/pause + prev/next
- Hover tooltip for framework atoms, guest atoms, and segmentation surface
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import trimesh
from ase import Atoms
from ase.data import covalent_radii
from ase.geometry import cellpar_to_cell
from ase.io import read
from skimage.measure import marching_cubes

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from ase.data import atomic_masses as atomic_masses_table

from mofography.data import read_raspa_pdb
import mofography as mgr
from mofography.helper import get_fractional_coordinates, interpolate_labels


CPK_COLORS = {
    "H": "#FFFFFF",
    "C": "#909090",
    "N": "#3050F8",
    "O": "#FF0D0D",
    "S": "#FFFF30",
    "P": "#FF8000",
    "F": "#90E050",
    "CL": "#1FF01F",
    "BR": "#A62929",
    "I": "#940094",
    "NA": "#AB5CF2",
    "K": "#8F40D4",
    "CA": "#3DFF00",
}


def _hex_to_rgb01(hex_color: str) -> list[float]:
    c = hex_color.lstrip("#")
    if len(c) != 6:
        return [0.8, 0.8, 0.8]
    return [int(c[0:2], 16) / 255.0, int(c[2:4], 16) / 255.0, int(c[4:6], 16) / 255.0]


def _serialize_movie(pdb_path: Path, max_frames: int | None) -> dict:
    movie = read_raspa_pdb(str(pdb_path))
    coords = movie["coords"]
    symbols = movie["symbols"]
    n_total = len(coords)
    n_frames = n_total if max_frames is None else min(n_total, max_frames)

    frames = []
    for i in range(n_frames):
        xyz = np.asarray(coords[i], dtype=float)
        syms = [str(s).upper() for s in symbols[i]]
        # Keep adsorbates visually distinct from framework: fixed pink spheres.
        colors = np.asarray([[1.0, 0.41, 0.71] for _ in syms], dtype=float)
        atomic_numbers = [Atoms(symbols=[s]).get_atomic_numbers()[0] for s in syms]
        radii = covalent_radii[np.asarray(atomic_numbers, dtype=int)] * 0.8
        frames.append(
            {
                "positions": xyz.reshape(-1).tolist(),
                "symbols": syms,
                "colors": colors.reshape(-1).tolist(),
                "radii": radii.astype(float).tolist(),
            }
        )

    return {
        "source": str(pdb_path),
        "n_frames": n_frames,
        "n_atoms": int(len(frames[0]["symbols"])) if frames else 0,
        "frames": frames,
    }


def _framework_vertex_hovertext(
    mesh: trimesh.Trimesh,
    ase_atoms,
    bond_rgba: tuple[int, int, int, int] = (127, 128, 129, 255),
) -> list[str]:
    """
    Match ``visualize.add_hoverdata_to_atoms_traces`` / Plotly Mesh3d hover lines
    (Id, Element, Mass, Position for atoms; bond indices for bond cylinders).
    """
    bonds = list(ase_atoms.info.get("covalent_bonds", []))
    n_atoms = int(ase_atoms.get_global_number_of_atoms())
    n_bonds = len(bonds)
    vc = np.asarray(mesh.visual.vertex_colors)
    if vc.size == 0 or vc.shape[0] != len(mesh.vertices):
        return [""] * len(mesh.vertices)
    if vc.shape[1] < 4:
        pad = np.full((vc.shape[0], 4 - vc.shape[1]), 255, dtype=vc.dtype)
        vc = np.hstack([vc, pad])
    vc_u8 = np.clip(np.rint(vc.astype(float) * 255.0) if vc.max() <= 1.0 else vc, 0, 255).astype(np.uint8)
    colors = ["rgba(" + ",".join(c.astype(str)) + ")" for c in vc_u8]
    vertexcolors = np.asarray(colors, dtype=object)
    bond_rgba_string = "rgba(" + ",".join(str(x) for x in bond_rgba) + ")"
    n_verts = len(mesh.vertices)
    sphere_mask = vertexcolors != bond_rgba_string
    bond_mask = vertexcolors == bond_rgba_string
    n_sphere_verts = int(np.sum(sphere_mask))
    n_bond_verts = int(np.sum(bond_mask))
    if n_atoms <= 0 or n_sphere_verts % n_atoms != 0:
        return [""] * n_verts
    vertices_per_sphere = n_sphere_verts // n_atoms
    if n_bonds <= 0:
        vertices_per_bond = 0
    elif n_bond_verts % n_bonds != 0:
        return [""] * n_verts
    else:
        vertices_per_bond = n_bond_verts // n_bonds

    atomic_numbers = ase_atoms.get_atomic_numbers()
    chemical_symbols = ase_atoms.get_chemical_symbols()
    masses = atomic_masses_table[atomic_numbers]
    positions = np.round(ase_atoms.get_positions(), decimals=4)
    pos_str = [" ".join(map(str, row)) for row in positions]

    vertex_atom_ids = np.repeat(np.arange(n_atoms), vertices_per_sphere)
    vertex_symbols = np.repeat(chemical_symbols, vertices_per_sphere)
    vertex_masses = np.repeat(masses, vertices_per_sphere)
    vertex_pos = np.repeat(pos_str, vertices_per_sphere)

    hover_atoms = [
        "Id: "
        + str(i)
        + "\nElement: "
        + str(s)
        + "\nMass: "
        + str(m)
        + "\nPosition: "
        + str(p)
        for i, s, m, p in zip(
            vertex_atom_ids, vertex_symbols, vertex_masses, vertex_pos, strict=False
        )
    ]
    bond_atom_ids = np.repeat(np.asarray(bonds, dtype=int), vertices_per_bond, axis=0)
    hover_bonds = [
        "Bond between atoms " + " and ".join(v.astype(str)) for v in bond_atom_ids
    ]
    return hover_atoms + hover_bonds


def _serialize_trimesh(mesh: trimesh.Trimesh) -> dict:
    v = np.asarray(mesh.vertices, dtype=float)
    f = np.asarray(mesh.faces, dtype=int)
    colors = None
    if hasattr(mesh, "visual") and hasattr(mesh.visual, "vertex_colors"):
        vc = np.asarray(mesh.visual.vertex_colors)
        if vc.ndim == 2 and vc.shape[0] == len(v):
            if vc.shape[1] >= 3:
                colors = vc[:, :3].astype(float) / 255.0
    if colors is None:
        colors = np.tile(np.asarray([[0.7, 0.7, 0.7]], dtype=float), (len(v), 1))
    return {
        "positions": v.reshape(-1).tolist(),
        "indices": f.reshape(-1).tolist(),
        "colors": colors.reshape(-1).tolist(),
        "n_vertices": int(len(v)),
        "n_triangles": int(len(f)),
    }


def _serialize_framework_mesh(cif_path: Path) -> dict:
    atoms = read(str(cif_path))
    atoms = mgr.update_with_covalent_bonds(atoms)
    mesh = mgr.mesh_for_atoms(
        atoms,
        kwargs_spheres=dict(scale=0.8),
        kwargs_cylinders=dict(radius=0.25),
    )
    payload = _serialize_trimesh(mesh)
    payload["vertex_hover"] = _framework_vertex_hovertext(mesh, atoms)
    if len(payload["vertex_hover"]) != payload["n_vertices"]:
        payload["vertex_hover"] = [""] * int(payload["n_vertices"])
    payload["source"] = str(cif_path)
    payload["n_atoms"] = int(len(atoms))
    return payload


def _serialize_cell_lines(cif_path: Path) -> dict:
    atoms = read(str(cif_path))
    A = cellpar_to_cell(atoms.cell.cellpar()).T
    a, b, c = A[:, 0], A[:, 1], A[:, 2]
    o = np.array([0.0, 0.0, 0.0])
    corners = [
        o,
        a,
        b,
        c,
        a + b,
        a + c,
        b + c,
        a + b + c,
    ]
    edges = [
        (0, 1), (0, 2), (0, 3),
        (1, 4), (2, 4),
        (2, 6), (3, 6),
        (3, 5), (1, 5),
        (4, 7), (5, 7), (6, 7),
    ]
    pts = []
    for i, j in edges:
        pts.extend([corners[i], corners[j]])
    arr = np.asarray(pts, dtype=float)
    return {"positions": arr.reshape(-1).tolist(), "n_segments": len(edges)}


def _serialize_segmentation_surface_mesh(
    cif_path: Path,
    spacing: float,
    isovalue: float,
    step_size: int,
    max_triangles: int,
) -> dict:
    atoms = read(str(cif_path))
    dgrid = mgr.dgrid_from_atoms_cpu(atoms, spacing=spacing, block_size=20000).astype(
        np.float32
    )
    seg = mgr.do_periodic_segmentation_rust(
        dgrid,
        atoms,
        spacing=spacing,
        mask_thickness=0.5,
        h=0.2,
        threshold_abs=1.5,
        min_distance=1.5,
        return_dgrid=True,
        return_padded_regions=False,
        extra_padding=5,
    )
    dgrid = seg["dgrid"]
    regions = seg["regions"]
    # Use one full isosurface mesh (smooth/continuous), then color vertices by
    # interpolated region labels to avoid region-boundary holes.
    verts, faces, _normals, _vals = marching_cubes(
        dgrid, level=isovalue, step_size=max(step_size, 1)
    )
    if len(faces) == 0:
        return {
            "positions": [],
            "indices": [],
            "colors": [],
            "labels": [],
            "n_vertices": 0,
            "n_triangles": 0,
        }
    cell = atoms.get_cell().array.T
    shape = np.asarray(dgrid.shape, dtype=float)
    verts_cart = np.dot(cell, (verts / shape).T).T
    mesh = trimesh.Trimesh(vertices=verts_cart, faces=faces, process=False)
    # Do NOT drop every Nth face: that creates visible holes/wireframe artifacts.
    # Keep mesh watertight; callers can raise step_size to reduce complexity.

    # Color vertices by region id with the same qualitative palette style.
    frac = get_fractional_coordinates(np.asarray(mesh.vertices, dtype=float), atoms)
    interp = interpolate_labels(regions)
    labels = interp(frac % 1).astype(int)
    cscale = [
        [102, 197, 204],
        [143, 200, 178],
        [184, 203, 152],
        [225, 206, 126],
        [246, 200, 113],
        [247, 171, 115],
        [240, 162, 152],
        [224, 173, 224],
    ]
    uniq = sorted([int(x) for x in np.unique(labels) if int(x) > 0])
    cmap = {rid: cscale[i % len(cscale)] for i, rid in enumerate(uniq)}
    colors = np.asarray(
        [cmap.get(int(r), [180, 180, 180]) for r in labels], dtype=float
    ) / 255.0
    return {
        "positions": np.asarray(mesh.vertices, dtype=float).reshape(-1).tolist(),
        "indices": np.asarray(mesh.faces, dtype=int).reshape(-1).tolist(),
        "colors": colors.reshape(-1).tolist(),
        "labels": labels.tolist(),
        "n_vertices": int(len(mesh.vertices)),
        "n_triangles": int(len(mesh.faces)),
    }


def _build_payload(
    pdb_path: Path,
    cif_path: Path,
    max_frames: int | None,
    spacing: float,
    isovalue: float,
    surf_step_size: int,
    max_triangles: int,
) -> dict:
    movie = _serialize_movie(pdb_path, max_frames=max_frames)
    framework = _serialize_framework_mesh(cif_path)
    cell = _serialize_cell_lines(cif_path)
    surface = _serialize_segmentation_surface_mesh(
        cif_path=cif_path,
        spacing=spacing,
        isovalue=isovalue,
        step_size=surf_step_size,
        max_triangles=max_triangles,
    )
    return {"movie": movie, "framework": framework, "surface": surface, "cell": cell}


def _html_template(payload_json: str, point_size: float, fps: float, eye, center, up) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Framework + Segmentation + MD Snapshots</title>
  <style>
    html, body {{
      margin: 0; padding: 0; width: 100%; height: 100%;
      overflow: hidden;
      background: radial-gradient(circle at center, #7a7a7a 0%, #111111 100%);
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: #efefef;
    }}
    #app {{ position: fixed; inset: 0; }}
    #hud {{
      position: fixed; left: 12px; right: 12px; bottom: 12px;
      background: rgba(20,20,20,0.75); border-radius: 10px; padding: 10px;
      display: grid; gap: 8px;
      backdrop-filter: blur(4px);
    }}
    #row {{ display: flex; align-items: center; gap: 8px; }}
    #frameSlider {{ flex: 1; }}
    .checkrow {{ display: flex; gap: 14px; align-items: center; flex-wrap: wrap; }}
    label {{ user-select: none; }}
    #tooltip {{
      position: fixed; pointer-events: none; display: none;
      background: rgba(0,0,0,0.85); border: 1px solid #666;
      border-radius: 8px; padding: 6px 8px; font-size: 12px;
      white-space: pre;
    }}
    button {{ background: #2e2e2e; color: #eee; border: 1px solid #555; border-radius: 6px; padding: 5px 9px; }}
  </style>
</head>
<body>
  <div id="app"></div>
  <div id="tooltip"></div>
  <div id="hud">
    <div id="row">
      <button id="playPause">Play</button>
      <button id="prevFrame">Prev</button>
      <button id="nextFrame">Next</button>
      <span id="status">Frame 0</span>
      <span id="meta"></span>
    </div>
    <div class="checkrow">
      <button id="pivotScene">Pivot: Scene center</button>
      <button id="pivotOrigin">Pivot: Origin (0,0,0)</button>
    </div>
    <div class="checkrow">
      <label><input type="checkbox" id="showCell" checked /> Unit cell</label>
      <label><input type="checkbox" id="showFramework" checked /> Framework mesh</label>
      <label><input type="checkbox" id="showSurface" checked /> Segmentation surface (regions)</label>
      <label><input type="checkbox" id="showGuests" checked /> Snapshot adsorbates (mesh spheres)</label>
    </div>
    <input id="frameSlider" type="range" min="0" max="0" step="1" value="0" />
  </div>

  <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
  <script>
    const DATA = {payload_json};
    const MOV = DATA.movie;
    const FW = DATA.framework;
    const SURF = DATA.surface;
    const CELL = DATA.cell;
    const app = document.getElementById("app");
    const slider = document.getElementById("frameSlider");
    const statusEl = document.getElementById("status");
    const metaEl = document.getElementById("meta");
    const playBtn = document.getElementById("playPause");
    const prevBtn = document.getElementById("prevFrame");
    const nextBtn = document.getElementById("nextFrame");
    const showFramework = document.getElementById("showFramework");
    const showCell = document.getElementById("showCell");
    const showSurface = document.getElementById("showSurface");
    const showGuests = document.getElementById("showGuests");
    const pivotSceneBtn = document.getElementById("pivotScene");
    const pivotOriginBtn = document.getElementById("pivotOrigin");
    const tooltip = document.getElementById("tooltip");

    const renderer = new THREE.WebGLRenderer({{ antialias: true, alpha: true }});
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    renderer.setSize(window.innerWidth, window.innerHeight);
    // Three.js r128 WebGLRenderer defaults (for side-by-side comparison with tuned lighting).
    renderer.outputEncoding = THREE.LinearEncoding;
    renderer.toneMapping = THREE.LinearToneMapping;
    renderer.toneMappingExposure = 1.0;
    app.appendChild(renderer.domElement);

    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.01, 2000);
    camera.position.set({eye[0]}, {eye[1]}, {eye[2]});
    camera.up.set({up[0]}, {up[1]}, {up[2]});
    camera.lookAt(new THREE.Vector3({center[0]}, {center[1]}, {center[2]}));
    const controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.target.set({center[0]}, {center[1]}, {center[2]});
    controls.update();

    // Soft ambient + key/fill directionals on opposite sides of the scene (parallel rays toward center).
    const ambient = new THREE.AmbientLight(0xffffff, 0.4);
    scene.add(ambient);
    const dirKey = new THREE.DirectionalLight(0xffffff, 0.85);
    dirKey.position.set(5, 10, 7.5);
    scene.add(dirKey);
    const dirFill = new THREE.DirectionalLight(0xffffff, 0.55);
    dirFill.position.set(-5, -9, -7.5);
    scene.add(dirFill);

    const frameworkGeo = new THREE.BufferGeometry();
    frameworkGeo.setAttribute("position", new THREE.Float32BufferAttribute(FW.positions, 3));
    frameworkGeo.setAttribute("color", new THREE.Float32BufferAttribute(FW.colors, 3));
    frameworkGeo.setIndex(FW.indices);
    frameworkGeo.computeVertexNormals();
    const frameworkMat = new THREE.MeshStandardMaterial({{
      vertexColors: true,
      transparent: true,
      opacity: 0.95,
      roughness: 0.38,
      metalness: 0.0,
      side: THREE.DoubleSide
    }});
    const frameworkMesh = new THREE.Mesh(frameworkGeo, frameworkMat);
    scene.add(frameworkMesh);

    const cellGeo = new THREE.BufferGeometry();
    cellGeo.setAttribute("position", new THREE.Float32BufferAttribute(CELL.positions, 3));
    const cellMat = new THREE.LineBasicMaterial({{ color: 0xffa500, linewidth: 2 }});
    const cellLines = new THREE.LineSegments(cellGeo, cellMat);
    scene.add(cellLines);

    const guestSphereGeo = new THREE.SphereBufferGeometry(1.0, 14, 14);
    const guestSphereMat = new THREE.MeshPhongMaterial({{
      color: 0xff69b4,
      shininess: 10,
      specular: 0x151515,
    }});
    const nGuest = Math.max(MOV.n_atoms, 1);
    const guestMesh = new THREE.InstancedMesh(guestSphereGeo, guestSphereMat, nGuest);
    guestMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    scene.add(guestMesh);

    const surfaceGeo = new THREE.BufferGeometry();
    surfaceGeo.setAttribute("position", new THREE.Float32BufferAttribute(SURF.positions, 3));
    surfaceGeo.setAttribute("color", new THREE.Float32BufferAttribute(SURF.colors, 3));
    surfaceGeo.setIndex(SURF.indices);
    surfaceGeo.computeVertexNormals();
    const surfaceMat = new THREE.MeshStandardMaterial({{
      vertexColors: true,
      transparent: true,
      opacity: 0.92,
      side: THREE.FrontSide,
      roughness: 0.34,
      metalness: 0.0
    }});
    const surfaceMesh = new THREE.Mesh(surfaceGeo, surfaceMat);
    scene.add(surfaceMesh);

    const sceneBox = new THREE.Box3();
    sceneBox.expandByObject(frameworkMesh);
    sceneBox.expandByObject(surfaceMesh);
    const sceneCenter = new THREE.Vector3();
    sceneBox.getCenter(sceneCenter);

    let frame = 0;
    let playing = false;
    let lastT = 0;
    const frameDtMs = 1000.0 / Math.max({fps}, 0.01);

    slider.max = String(Math.max(MOV.n_frames - 1, 0));
    metaEl.textContent = `Framework atoms: ${{FW.n_atoms}} | Snapshot atoms: ${{MOV.n_atoms}} | Frames: ${{MOV.n_frames}} | Surface tris: ${{SURF.n_triangles}}`;
    controls.target.copy(sceneCenter);
    controls.update();

    function setFrame(i) {{
      frame = Math.max(0, Math.min(i, MOV.n_frames - 1));
      const fr = MOV.frames[frame];
      const dummy = new THREE.Object3D();
      for (let a = 0; a < MOV.n_atoms; a++) {{
        const x = fr.positions[a * 3], y = fr.positions[a * 3 + 1], z = fr.positions[a * 3 + 2];
        const r = fr.radii[a] || 0.6;
        dummy.position.set(x, y, z);
        dummy.scale.set(r, r, r);
        dummy.updateMatrix();
        guestMesh.setMatrixAt(a, dummy.matrix);
      }}
      guestMesh.count = MOV.n_atoms;
      guestMesh.instanceMatrix.needsUpdate = true;
      statusEl.textContent = `Frame ${{frame}}`;
      slider.value = String(frame);
    }}

    slider.addEventListener("input", () => setFrame(parseInt(slider.value, 10)));
    prevBtn.addEventListener("click", () => setFrame(frame - 1));
    nextBtn.addEventListener("click", () => setFrame(frame + 1));
    playBtn.addEventListener("click", () => {{
      playing = !playing;
      playBtn.textContent = playing ? "Pause" : "Play";
    }});
    showCell.addEventListener("change", () => cellLines.visible = showCell.checked);
    showFramework.addEventListener("change", () => frameworkMesh.visible = showFramework.checked);
    showSurface.addEventListener("change", () => surfaceMesh.visible = showSurface.checked);
    showGuests.addEventListener("change", () => guestMesh.visible = showGuests.checked);
    pivotSceneBtn.addEventListener("click", () => {{
      controls.target.copy(sceneCenter);
      controls.update();
    }});
    pivotOriginBtn.addEventListener("click", () => {{
      controls.target.set(0, 0, 0);
      controls.update();
    }});

    const raycaster = new THREE.Raycaster();
    raycaster.params.Points.threshold = 0.2;
    const mouse = new THREE.Vector2();
    function onPointerMove(ev) {{
      const rect = renderer.domElement.getBoundingClientRect();
      mouse.x = ((ev.clientX - rect.left) / rect.width) * 2 - 1;
      mouse.y = -((ev.clientY - rect.top) / rect.height) * 2 + 1;
      raycaster.setFromCamera(mouse, camera);
      const fwHits = frameworkMesh.visible ? raycaster.intersectObject(frameworkMesh) : [];
      const gHits = guestMesh.visible ? raycaster.intersectObject(guestMesh) : [];
      const sHits = surfaceMesh.visible ? raycaster.intersectObject(surfaceMesh) : [];
      const all = [...fwHits, ...gHits, ...sHits].sort((a, b) => a.distance - b.distance);
      if (!all.length) {{
        tooltip.style.display = "none";
        return;
      }}
      const hit = all[0];
      let txt = "";
      if (hit.object === guestMesh) {{
        const idx = hit.instanceId ?? -1;
        const fr = MOV.frames[frame];
        const p = fr.positions.slice(idx * 3, idx * 3 + 3);
        const s = fr.symbols[idx] || "?";
        txt = `trace: adsorbate spheres\\nframe: ${{frame}}\\natom: ${{idx}}\\nelement: ${{s}}\\nxyz: (${{p[0].toFixed(3)}}, ${{p[1].toFixed(3)}}, ${{p[2].toFixed(3)}})`;
      }} else if (hit.object === frameworkMesh) {{
        const vh = FW.vertex_hover;
        let hoverLine = "";
        let vi = -1;
        if (hit.face) {{
          vi = hit.face.a;
        }} else if (hit.faceIndex != null && FW.indices && hit.faceIndex >= 0) {{
          vi = FW.indices[hit.faceIndex * 3] ?? -1;
        }}
        if (vh && Array.isArray(vh) && vi >= 0 && vi < vh.length && vh[vi]) {{
          hoverLine = String(vh[vi]);
        }}
        if (hoverLine) {{
          txt = `trace: framework mesh\\n${{hoverLine}}`;
        }} else {{
          txt = `trace: framework mesh\\nxyz: (${{hit.point.x.toFixed(3)}}, ${{hit.point.y.toFixed(3)}}, ${{hit.point.z.toFixed(3)}})`;
        }}
      }} else if (hit.object === surfaceMesh) {{
        const fi = hit.faceIndex ?? -1;
        let region = "?";
        if (fi >= 0) {{
          const ia = SURF.indices[fi * 3] ?? -1;
          if (ia >= 0) region = String(SURF.labels[ia] ?? "?");
        }}
        txt = `trace: segmentation surface\\nregion label: ${{region}}\\nxyz: (${{hit.point.x.toFixed(3)}}, ${{hit.point.y.toFixed(3)}}, ${{hit.point.z.toFixed(3)}})`;
      }}
      tooltip.style.display = "block";
      tooltip.style.left = `${{ev.clientX + 12}}px`;
      tooltip.style.top = `${{ev.clientY + 12}}px`;
      tooltip.textContent = txt;
    }}
    renderer.domElement.addEventListener("pointermove", onPointerMove);
    renderer.domElement.addEventListener("pointerleave", () => tooltip.style.display = "none");

    function animate(t) {{
      requestAnimationFrame(animate);
      if (playing && MOV.n_frames > 1 && (t - lastT) >= frameDtMs) {{
        setFrame((frame + 1) % MOV.n_frames);
        lastT = t;
      }}
      controls.update();
      renderer.render(scene, camera);
    }}

    window.addEventListener("resize", () => {{
      camera.aspect = window.innerWidth / window.innerHeight;
      camera.updateProjectionMatrix();
      renderer.setSize(window.innerWidth, window.innerHeight);
    }});

    cellLines.visible = showCell.checked;
    frameworkMesh.visible = showFramework.checked;
    surfaceMesh.visible = showSurface.checked;
    guestMesh.visible = showGuests.checked;
    setFrame(0);
    requestAnimationFrame(animate);
  </script>
</body>
</html>
"""


def main() -> None:
    p = argparse.ArgumentParser(
        description="Generate interactive Three.js HTML for framework + segmentation + MD snapshots."
    )
    p.add_argument(
        "--pdb",
        type=Path,
        default=Path("mofography/examples/Movie_LTA96_1.1.1_298.150000_0.000000_component_K_0.pdb"),
        help="Path to RASPA/MD PDB movie snapshots.",
    )
    p.add_argument(
        "--cif",
        type=Path,
        default=Path("examples/TestCIFs/LTA96.cif"),
        help="Framework CIF used for framework trace + segmentation surface.",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("outputs/md_pdb_threejs_framework_surface_viewer.html"),
        help="Output HTML path.",
    )
    p.add_argument("--max-frames", type=int, default=200, help="Optional frame cap.")
    p.add_argument("--point-size", type=float, default=0.16, help="Three.js point size.")
    p.add_argument("--fps", type=float, default=20.0, help="Playback FPS in viewer.")
    p.add_argument("--spacing", type=float, default=0.25, help="Segmentation dgrid spacing.")
    p.add_argument("--surface-isovalue", type=float, default=1.0, help="Surface isovalue.")
    p.add_argument(
        "--surface-step-size",
        type=int,
        default=2,
        help="Marching cubes step size (higher=faster, coarser).",
    )
    p.add_argument(
        "--max-triangles",
        type=int,
        default=120000,
        help="Triangle cap for exported surface mesh.",
    )
    p.add_argument("--eye", nargs=3, type=float, default=(2.0, 0.0, 1.35))
    p.add_argument("--center", nargs=3, type=float, default=(0.0, 0.0, 0.0))
    p.add_argument("--up", nargs=3, type=float, default=(0.0, 0.0, 1.0))
    args = p.parse_args()

    pdb = args.pdb.expanduser().resolve()
    if not pdb.is_file():
        raise SystemExit(f"PDB not found: {pdb}")
    cif = args.cif.expanduser().resolve()
    if not cif.is_file():
        raise SystemExit(f"CIF not found: {cif}")

    payload = _build_payload(
        pdb_path=pdb,
        cif_path=cif,
        max_frames=args.max_frames,
        spacing=float(args.spacing),
        isovalue=float(args.surface_isovalue),
        surf_step_size=int(args.surface_step_size),
        max_triangles=int(args.max_triangles),
    )
    html = _html_template(
        payload_json=json.dumps(payload),
        point_size=float(args.point_size),
        fps=float(args.fps),
        eye=tuple(args.eye),
        center=tuple(args.center),
        up=tuple(args.up),
    )

    out = args.out.expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
