#!/usr/bin/env python3
"""
NU-1000 Periodic Segmentation Visualization (Simplified)

Focuses on pore network visualization without complex framework traces
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
    
    # Compute distance grid
    print(f"\n📏 Computing distance grid...")
    dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=SPACING)
    dgrid = mgr.zoom_grid(dgrid, zoom_factor=1).compute()
    dgrid = mgr.periodic_gaussian_filter(dgrid, sigma=SPACING*2)
    dgrid = dgrid.astype(np.float32)
    
    print(f"✓ Distance grid shape: {dgrid.shape}")
    
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
    connections = results['connections']
    all_maxima = results['all_maxima']
    unique_regions = results['unique_regions']
    pbc_groups = results['pbc_groups']
    
    print(f"✓ Found {len(maxima_pbc)} unique periodic maxima")
    print(f"  Found {len(unique_regions)} unique regions")
    print(f"  Found {len(pbc_groups)} PBC groups")
    print(f"  Found {len(connections)} region connections")
    
    # Trim connections
    print(f"\n✂️  Trimming connections...")
    connections_trimmed = np.array([
        c for c in connections
        if np.any(np.isin(c[0], unique_regions))
    ])
    print(f"  Original: {len(connections)} → Trimmed: {len(connections_trimmed)}")
    
    # Create pore graph
    print(f"\n🕸️  Creating pore graph...")
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
        print(f"  Warning: {e}")
        lcd, pld, dfsd = 0.0, 0.0, 0.0
    
    # Create simplified visualization
    print(f"\n🎨 Creating 3D visualization...")
    
    cscale = px.colors.qualitative.Pastel
    traces = []
    
    # 1. Unit cell (wireframe box)
    print(f"  ✓ Adding unit cell...")
    cell = mat_atoms.get_cell()
    a, b, c = cell
    
    # Create wireframe box
    def create_box_edges(a, b, c):
        edges = []
        # Define 8 corners
        corners = np.array([
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]
        ])
        # Transform corners
        corners = corners @ [a, b, c]
        
        # Define 12 edges
        edge_pairs = [
            (0,1), (1,2), (2,3), (3,0),  # bottom face
            (4,5), (5,6), (6,7), (7,4),  # top face
            (0,4), (1,5), (2,6), (3,7)   # vertical edges
        ]
        
        for i, j in edge_pairs:
            edges.extend([
                corners[i][0], corners[j][0], None,
                corners[i][1], corners[j][1], None,
                corners[i][2], corners[j][2], None
            ])
        return edges
    
    box_edges = create_box_edges(a, b, c)
    
    cell_trace = go.Scatter3d(
        x=box_edges[0::3],
        y=box_edges[1::3],
        z=box_edges[2::3],
        mode='lines',
        line=dict(color='white', width=2),
        name='Unit Cell',
        showlegend=True
    )
    traces.append(cell_trace)
    
    # 2. Segmented surface with PBC (isosurface)
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
    
    # 3. Local maxima with PBC (pore centers)
    print(f"  ✓ Adding pore centers (maxima PBC)...")
    # Transform maxima to Cartesian coordinates
    maxima_cartesian = maxima_pbc @ [a, b, c]
    
    maxima_trace = go.Scatter3d(
        x=maxima_cartesian[:, 0],
        y=maxima_cartesian[:, 1],
        z=maxima_cartesian[:, 2],
        mode='markers',
        marker=dict(
            size=8,
            color=list(range(len(maxima_pbc))),
            colorscale=cscale,
            opacity=0.9,
            line=dict(color='white', width=1)
        ),
        text=[f'Pore {i+1}' for i in range(len(maxima_pbc))],
        name='Pore centers (PBC)',
        showlegend=True
    )
    traces.append(maxima_trace)
    
    # 4. Pore graph/connections (lines between connected pores)
    print(f"  ✓ Adding pore graph connections...")
    for i, connection in enumerate(connections_trimmed[:50]):  # Limit to 50 for performance
        region_pair = connection[0]
        if len(region_pair) >= 2:
            # Get the maxima indices for these regions
            region1_idx, region2_idx = region_pair[0], region_pair[1]
            
            # Find which maxima correspond to these regions
            max1_idx = region1_idx - 1
            max2_idx = region2_idx - 1
            
            if 0 <= max1_idx < len(all_maxima) and 0 <= max2_idx < len(all_maxima):
                max1_pos = all_maxima[max1_idx] @ [a, b, c]
                max2_pos = all_maxima[max2_idx] @ [a, b, c]
                
                conn_trace = go.Scatter3d(
                    x=[max1_pos[0], max2_pos[0], None],
                    y=[max1_pos[1], max2_pos[1], None],
                    z=[max1_pos[2], max2_pos[2], None],
                    mode='lines',
                    line=dict(color='rgba(255,255,255,0.3)', width=1),
                    showlegend=False,
                    hoverinfo='skip'
                )
                traces.append(conn_trace)
    
    # Create figure
    print(f"\n🖼️  Creating figure...")
    fig = go.Figure(data=traces)
    
    fig.update_layout(
        scene=dict(
            xaxis=dict(title='X (Å)', backgroundcolor='rgb(24, 24, 27)'),
            yaxis=dict(title='Y (Å)', backgroundcolor='rgb(24, 24, 27)'),
            zaxis=dict(title='Z (Å)', backgroundcolor='rgb(24, 24, 27)'),
            aspectmode='data',
            bgcolor='rgb(24, 24, 27)'
        ),
        title=dict(
            text=f"<b>NU-1000: Pore Network Analysis</b><br>" +
                  f"<sub>LCD: {lcd:.2f} Å | PLD: {pld:.2f} Å | DFSD: {dfsd:.2f} Å | " +
                  f"Pores: {len(maxima_pbc)} | Connections: {len(connections_trimmed)}</sub>",
            x=0.5,
            xanchor='center',
            font=dict(size=16, color='white')
        ),
        paper_bgcolor='rgb(24, 24, 27)',
        plot_bgcolor='rgb(24, 24, 27)',
        font=dict(color='white'),
        showlegend=True,
        legend=dict(
            xanchor="left",
            yanchor="top",
            x=0.02,
            y=0.98,
            bgcolor='rgba(0,0,0,0.5)',
            bordercolor='white',
            borderwidth=1
        ),
        margin=dict(l=0, r=0, b=0, t=80)
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
        f.write(f"Unit cell dimensions: {mat_atoms.get_cell_lengths_and_angles()[:3]}\n\n")
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
        f.write(f"\nVisualization Features:\n")
        f.write(f"  ✓ Unit cell (wireframe box)\n")
        f.write(f"  ✓ Segmented surface with PBC (pore surfaces)\n")
        f.write(f"  ✓ Local maxima with PBC (pore centers)\n")
        f.write(f"  ✓ Pore graph/connections (lines between pores)\n")
    
    print(f"✓ Saved: nu1000_analysis_summary.txt")
    
    print(f"\n" + "=" * 60)
    print("✅ ANALYSIS COMPLETE!")
    print("=" * 60)
    print(f"\n📊 Results:")
    print(f"  HTML visualization: {output_html.absolute()}")
    print(f"  Summary: nu1000_analysis_summary.txt")
    print(f"\n🌐 Open in browser:")
    print(f"  open {output_html.absolute()}")
    print(f"\n" + "=" * 60)
    
    return output_html.absolute()

if __name__ == "__main__":
    main()