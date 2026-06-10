#!/usr/bin/env python3
"""
NU-1000 Periodic Segmentation Visualization

Complete visualization of NU-1000 structure showing:
1. Unit cell
2. Segmented surface with PBC
3. Pore graph/connections
4. Local maxima with PBC
"""

import numpy as np
from pathlib import Path
import mofography as mgr
from ase.io import read
import plotly.graph_objects as go
import plotly.express as px

# Parameters
SPACING = 0.5
MASK_THICKNESS = 1.0
H_MAXIMA = 0.2
THRESHOLD_ABS = 1.5
MINIMUM_WINDOW_SEPARATION = 1.5
WALL_THICKNESS = 0.08
EXTRA_PADDING = 10
ISOVALUE = 1.5

def main():
    print("=" * 60)
    print("NU-1000 Periodic Segmentation Analysis")
    print("=" * 60)
    
    # Load structure
    print("\n📁 Loading NU-1000 structure...")
    cif_path = Path("mofography/examples/TestCIFs/NU-1000_S1.cif")
    mat_atoms = read(str(cif_path))
    
    print(f"✓ Loaded: {mat_atoms.get_chemical_formula()}")
    print(f"  Atoms: {len(mat_atoms)}")
    print(f"  Cell: {mat_atoms.get_cell_lengths_and_angles()[:3]}")
    
    # Compute distance grid
    print(f"\n📏 Computing distance grid (spacing={SPACING}Å)...")
    dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=SPACING)
    dgrid = mgr.zoom_grid(dgrid, zoom_factor=1).compute()
    dgrid = mgr.periodic_gaussian_filter(dgrid, sigma=SPACING*2)
    dgrid = dgrid.astype(np.float32)  # Critical for Rust backend
    
    print(f"✓ Distance grid shape: {dgrid.shape}")
    print(f"  Range: {dgrid.min():.2f} to {dgrid.max():.2f} Å")
    
    # Run periodic segmentation
    print(f"\n🎯 Running periodic segmentation...")
    results = mgr.do_periodic_segmentation_rust(
        dgrid,
        mat_atoms,
        spacing=SPACING,
        mask_thickness=MASK_THICKNESS,
        h=H_MAXIMA,
        threshold_abs=THRESHOLD_ABS,
        minimum_window_separation=MINIMUM_WINDOW_SEPARATION,
        wall_thickness=WALL_THICKNESS,
        extra_padding=EXTRA_PADDING,
        return_padded_regions=True,
        compute_connections=True
    )
    
    # Extract results
    maxima_pbc = results['maxima_pbc']
    regions_pbc = results['regions_pbc']
    regions = results['regions']
    connections = results['connections']
    all_maxima = results['all_maxima']
    unique_regions = results['unique_regions']
    pbc_groups = results['pbc_groups']
    
    print(f"✓ Found {len(maxima_pbc)} unique periodic maxima")
    print(f"  Found {len(unique_regions)} unique regions")
    print(f"  Found {len(pbc_groups)} PBC groups")
    print(f"  Found {len(connections)} region connections")
    
    # Trim connections
    print(f"\n✂️  Trimming connections to unit cell...")
    connections_trimmed = np.array([
        c for c in connections
        if np.any(np.isin(c[0], unique_regions))
    ])
    print(f"  Original: {len(connections)} → Trimmed: {len(connections_trimmed)}")
    
    # Create pore graph
    print(f"\n🕸️  Creating pore graph (RAG)...")
    rag = mgr.rag_from_connections(connections_trimmed, all_maxima)
    print(f"  Nodes: {rag.number_of_nodes()}")
    print(f"  Edges: {rag.number_of_edges()}")
    
    # Compute pore diameters
    print(f"\n📊 Computing pore diameters...")
    try:
        lcd, pld, dfsd, _ = mgr.get_pore_dia_mgr(
            mat_atoms, dgrid=dgrid, 
            kwargs=dict(spacing=SPACING)
        )
        print(f"  LCD: {lcd:.3f} Å")
        print(f"  PLD: {pld:.3f} Å")
        print(f"  DFSD: {dfsd:.3f} Å")
    except Exception as e:
        print(f"  Warning: Pore diameter calculation failed: {e}")
        lcd, pld, dfsd = 0.0, 0.0, 0.0
    
    # Create visualization
    print(f"\n🎨 Creating 3D visualization...")
    
    # Color scale
    cscale = px.colors.qualitative.Pastel
    
    # Update structure with covalent bonds
    mat_updated = mgr.update_with_covalent_bonds(mat_atoms)
    
    # Create traces
    traces = []
    
    # 1. Unit cell
    print(f"  ✓ Adding unit cell...")
    cell_trace = mgr.plot_cell(mat_atoms)
    traces.append(cell_trace)
    
    # 2. Framework
    print(f"  ✓ Adding framework...")
    frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)
    traces.append(frame_trace)
    
    # 3. Segmented surface with PBC
    print(f"  ✓ Adding segmented surface (PBC merged)...")
    surface_trace = mgr.plot_surface_dgrid(
        dgrid,
        isovalue=ISOVALUE,
        region_labels=regions_pbc,
        ase_atoms=mat_atoms,
        cscale=cscale,
        name='Pores (PBC merged)',
        stylize=True,
        opacity=0.8
    )
    traces.append(surface_trace)
    
    # 4. Local maxima with PBC
    print(f"  ✓ Adding pore centers (maxima PBC)...")
    maxima_trace = mgr.plot_maxima(
        dgrid, maxima_pbc, mat_atoms,
        name='Pore centers (PBC)',
        show_spheres=False,
        marker_colorscale=cscale
    )
    traces.append(maxima_trace)
    
    # 5. Pore graph/connections
    print(f"  ✓ Adding pore graph connections...")
    connections_traces = mgr.plot_connections_as_it_is(
        connections_trimmed,
        all_maxima,
        dgrid.shape,
        mat_atoms,
        cscale=cscale,
        plot_disconnected=False,
        show_spheres=False
    )
    traces.extend(connections_traces)
    
    # Create figure
    print(f"\n🖼️  Creating figure...")
    fig = go.Figure(
        data=traces,
        layout=mgr.dark_layout
    )
    
    fig.update_layout(
        scene_aspectmode='data',
        title=dict(
            text=f"<b>NU-1000: Pore Network Analysis</b><br>" +
                  f"<sub>LCD: {lcd:.2f} Å | PLD: {pld:.2f} Å | DFSD: {dfsd:.2f} Å | " +
                  f"Pores: {len(maxima_pbc)} | Connections: {len(connections_trimmed)}</sub>",
            x=0.5,
            xanchor='center',
            font=dict(size=16, color='white')
        ),
        scene=dict(
            xaxis=dict(title='X (Å)'),
            yaxis=dict(title='Y (Å)'),
            zaxis=dict(title='Z (Å)')
        ),
        paper_bgcolor='rgb(24, 24, 27)',
        plot_bgcolor='rgb(24, 24, 27)'
    )
    
    # Save HTML
    output_html = Path("nu1000_periodic_segmentation.html")
    print(f"\n💾 Saving to {output_html}...")
    fig.write_html(str(output_html), include_plotlyjs='cdn')
    print(f"✓ Saved: {output_html.absolute()}")
    
    # Save summary
    print(f"\n📋 Saving analysis summary...")
    with open("nu1000_analysis_summary.txt", "w") as f:
        f.write("NU-1000 Periodic Segmentation Analysis\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Structure: {mat_atoms.get_chemical_formula()}\n")
        f.write(f"Atoms: {len(mat_atoms)}\n")
        f.write(f"Unit cell: {mat_atoms.get_cell_lengths_and_angles()[:3]}\n\n")
        f.write("Parameters:\n")
        f.write(f"  Grid spacing: {SPACING} Å\n")
        f.write(f"  Mask thickness: {MASK_THICKNESS} Å\n")
        f.write(f"  H-maxima: {H_MAXIMA}\n")
        f.write(f"  Threshold abs: {THRESHOLD_ABS} Å\n\n")
        f.write("Results:\n")
        f.write(f"  Unique periodic maxima: {len(maxima_pbc)}\n")
        f.write(f"  Unique regions: {len(unique_regions)}\n")
        f.write(f"  PBC groups: {len(pbc_groups)}\n")
        f.write(f"  Total connections: {len(connections)}\n")
        f.write(f"  Trimmed connections: {len(connections_trimmed)}\n")
        f.write(f"  LCD: {lcd:.3f} Å\n")
        f.write(f"  PLD: {pld:.3f} Å\n")
        f.write(f"  DFSD: {dfsd:.3f} Å\n")
    
    print(f"✓ Saved: nu1000_analysis_summary.txt")
    
    # Save NPZ
    print(f"\n💾 Saving NPZ results...")
    np.savez_compressed(
        "nu1000_segmentation_results.npz",
        dgrid=dgrid,
        maxima_pbc=maxima_pbc,
        regions_pbc=regions_pbc,
        regions=regions,
        connections=connections_trimmed,
        all_maxima=all_maxima,
        unique_regions=unique_regions,
        pbc_groups=pbc_groups
    )
    print(f"✓ Saved: nu1000_segmentation_results.npz")
    
    print(f"\n" + "=" * 60)
    print("✅ ANALYSIS COMPLETE!")
    print("=" * 60)
    print(f"\n📊 Results:")
    print(f"  HTML visualization: {output_html.absolute()}")
    print(f"  Summary: nu1000_analysis_summary.txt")
    print(f"  Data: nu1000_segmentation_results.npz")
    print(f"\n🌐 Open the HTML file in your browser:")
    print(f"  open {output_html.absolute()}")
    print(f"\n" + "=" * 60)
    
    return output_html.absolute()

if __name__ == "__main__":
    main()