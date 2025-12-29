# Tutorial Creation Guidelines for MOFography

> **Purpose**: This document captures best practices and lessons learned from creating the Periodic Segmentation Quick Start tutorial. Use this as a reference when creating new tutorials or updating existing ones.

---

## Table of Contents

1. [Philosophy](#philosophy)
2. [Format Decision: Markdown vs Jupyter](#format-decision-markdown-vs-jupyter)
3. [Key Tutorial Sections](#key-tutorial-sections)
4. [Writing Style Guidelines](#writing-style-guidelines)
5. [Code Patterns](#code-patterns)
6. [Image Generation](#image-generation)
7. [Common Pitfalls](#common-pitfalls)
8. [Validation Checklist](#validation-checklist)

---

## Philosophy

### What Makes a Good Tutorial?

A good MOFography tutorial should be:

1. **Beginner-Friendly** - Assumes no prior knowledge of the package
2. **Self-Contained** - Can be followed without external references
3. **Engaging** - Uses conversational tone, emojis, and clear explanations
4. **Visual** - Includes images that actually work and show expected outputs
5. **Honest About Scope** - Clearly states what it covers and what it doesn't
6. **Practical** - Focuses on getting started quickly, not exhaustive documentation

### The "Why Before How" Principle

Before explaining HOW to do something, always explain WHY it matters:
- What problem does this solve?
- Why is this approach better than alternatives?
- What insights does this enable?

---

## Format Decision: Markdown vs Jupyter

### Recommendation: **Use Markdown for Tutorials**

After extensive testing, we found markdown files are superior to Jupyter notebooks for tutorials:

**Advantages of Markdown:**
- ✅ Better readability - plain text, no UI clutter
- ✅ Version control friendly - clean diffs
- ✅ Easy to edit - any text editor works
- ✅ Faster navigation - scroll through naturally
- ✅ Can still include Python code blocks with syntax highlighting
- ✅ Images link reliably
- ✅ Can be rendered by GitHub/GitLab natively

**When to Use Jupyter:**
- Interactive notebooks for data exploration
- Complex multi-step analyses
- When you need to actually execute code cells

**When to Use Markdown:**
- Documentation tutorials
- Quick start guides
- Reference material
- Anything meant to be read, not executed

---

## Key Tutorial Sections

Based on the Periodic Segmentation Quick Start tutorial, here's the template:

### 1. Header with Hooks

```markdown
# 🌟 [Engaging Title]

> **[One-sentence summary of what the tutorial covers]**

---

## 📚 Citation
```

**Include:**
- Catchy title with emoji
- Subtitle explaining the purpose
- Paper citation (if applicable) with bibtex

### 2. Conceptual Introduction

```markdown
## 🎯 What is [Topic]?

[Explain the core concept in simple terms]

### Why is this Powerful? 🚀

[3-5 bullet points explaining benefits]
```

**Best Practices:**
- Start simple, then add depth
- Use analogies where helpful
- Focus on insights enabled, not just features
- Use emojis to break up text

### 3. Tutorial Scope

```markdown
## ⚠️ Tutorial Scope

**This is a [SCOPE] tutorial**, not a comprehensive guide!

✅ **What this tutorial covers:**
- [Key topics covered]

❌ **What this tutorial does NOT cover:**
- [Topics intentionally omitted]

🔗 **For deeper understanding**, see [links to related tutorials].
```

**Why This Matters:**
- Manages user expectations
- Prevents tutorial bloat
- Guides users to appropriate resources

### 4. Context Setting

```markdown
## 🔄 The Old Workflow (Legacy)

[Explain the old approach and its problems]

**Problems with the old approach:**
- ❌ [Problem 1]
- ❌ [Problem 2]

**The New Way (Current):**
- ✅ [Solution 1]
- ✅ [Solution 2]
```

**When to Include:**
- When introducing a new workflow that replaces an old one
- When comparing approaches helps understanding

### 5. Motivation

```markdown
## 🌟 Why Do We Need [Feature]?

[Explain the problem in relatable terms]

![Descriptive image](./path/to/image.png)

*Caption explaining what the image shows*
```

**Best Practices:**
- Use concrete examples
- Show, don't just tell (use images!)
- Keep captions concise but informative

### 6. Learning Objectives

```markdown
## 🎓 What You'll Learn

By the end of this tutorial, you'll be able to:

1. ✅ [Skill 1]
2. ✅ [Skill 2]
3. ✅ [Skill 3]
```

### 7. Prerequisites

```markdown
## 🔧 Prerequisites

- [Requirement 1]
- [Requirement 2]
- (Optional) [Requirement 3]
```

### 8. The Workflow (Step-by-Step)

```markdown
## 📖 The Workflow in a Nutshell

```
[Visual representation of workflow steps]
```

---

## Step 1: [Title] [🔧/📦/🏗️/etc.]

[Explanation of what this step does and why]

```python
[code here]
```

**What is this?**
- [Brief explanation]

---

## Step 2: [Title] [🔧/📦/🏗️/etc.]

[...]
```

**Best Practices:**
- Use emojis in step titles for visual scanning
- Provide code that can be copy-pasted
- Follow code with "What just happened?" or "What is this?" sections
- Keep steps focused on one concept each

### 9. Parameter Guide

When introducing parameters, use a table:

```markdown
**Parameter Guide:**

| Parameter | What it does | Typical range |
|-----------|--------------|---------------|
| `param1`  | Description  | 0.1 - 1.0     |
| `param2`  | Description  | 1 - 10        |

**Tip:** Start with defaults and adjust only if results don't look right!
```

### 10. Output Explanation

```markdown
**The Output Arrays:**

| Output | What is it? | Why care? |
|--------|-------------|-----------|
| `var1` | Description | Importance |
| `var2` | Description | Importance |
```

### 11. Troubleshooting

```markdown
## 🛠️ Troubleshooting

| Issue | Solution |
|-------|----------|
| **Problem** | Fix |
| **Problem** | Fix |
```

### 12. Next Steps

```markdown
## 🚀 Next Steps

### Learn More in Detail:
- [Link 1] - Description
- [Link 2] - Description

### Explore Basic Tutorials:
- `notebook1.ipynb` - Description
- `notebook2.ipynb` - Description

### Explore Advanced Tutorials:
- `notebook3.ipynb` - Description
```

### 13. Summary

```markdown
## 📝 Summary

| Step | Function | Output |
|------|----------|--------|
| 1. [Step] | [Function] | [Output] |
| 2. [Step] | [Function] | [Output] |

**Key Takeaway:** [One-sentence summary of the main point]
```

### 14. Bonus Section

```markdown
## 🎁 Bonus: [Title]

[Additional helpful content]
```

---

## Writing Style Guidelines

### Tone and Voice

- **Conversational but professional** - Like a knowledgeable colleague explaining something
- **Enthusiastic** - Show why this is exciting/interesting
- **Direct** - Say what you mean, avoid wordiness
- **Second-person** - Use "you" to address the reader directly

### Emojis Usage

Use emojis to:
- Break up walls of text
- Provide visual anchors for scanning
- Indicate section types (warnings, tips, code, etc.)

**Common emoji patterns:**
- 🌟 For main section titles
- ✅ For checklists and positive outcomes
- ❌ For problems/warnings
- ⚠️ For warnings/cautions
- 🔧 For setup/configuration
- 📦 For imports/packages
- 🏗️ For structure/objects
- 📏 For measurements/grids
- 🎯 For objectives/targets
- 🚀 For next steps/advanced topics
- 🛠️ For troubleshooting
- 📝 For summaries
- 🎁 For bonus content
- 💡 For tips

**Don't overuse it** - Not every paragraph needs an emoji

### Formatting Best Practices

1. **Code blocks**: Always specify language for syntax highlighting
   ```python
   code here
   ```

2. **Inline code**: Use backticks for variable/function names
   - Use `do_periodic_segmentation_rust()` for periodic segmentation
   - The `dgrid` variable contains the distance grid

3. **Bold for emphasis**: Use sparingly
   - **CRITICAL**: Convert to float32
   - **Key insight**: Pores don't respect boundaries

4. **Tables**: Use for structured information
   - Parameter guides
   - Output descriptions
   - Troubleshooting tips

---

## Code Patterns

### Path Handling

**Use pathlib instead of os.path:**

```python
from pathlib import Path

# Get script directory
script_path = Path(__file__).resolve()
script_dir = script_path.parent

# Find repo root
if script_dir.name == 'docs':
    repo_root = script_dir.parent
else:
    repo_root = script_dir

# Build paths with proper handling
cif_name = 'P1-IRMOF-1.cif'
possible_cif_paths = [
    repo_root / 'examples' / 'TestCIFs' / cif_name,
    repo_root / 'mofography' / 'examples' / 'TestCIFs' / cif_name,
]

# Always use expanduser() for tilde paths
path_to_cif = None
for cif_path in possible_cif_paths:
    expanded_path = cif_path.expanduser()
    if expanded_path.exists():
        path_to_cif = expanded_path
        break
```

**Why:**
- More readable
- Handles cross-platform paths better
- `expanduser()` handles `~` symbol correctly

### Environment Variables

```python
# Allow user to override with environment variable
if 'MOFOGRAPHY_TEST_CIF' in os.environ:
    path_to_cif = Path(os.environ['MOFOGRAPHY_TEST_CIF']).expanduser()
```

### Critical Type Conversions

**ALWAYS document dtype requirements:**

```python
# Compute distance grid
dgrid = mgr.dgrid_from_atoms_cpu(mat_atoms, spacing=spacing, block_size=50000)
dgrid = mgr.zoom_grid(dgrid, zoom_factor=1).compute()
dgrid = mgr.periodic_gaussian_filter(dgrid, sigma=spacing*2)

# CRITICAL: Convert to float32 for Rust backend
original_dtype = dgrid.dtype
dgrid = dgrid.astype(np.float32)

print(f"  Original dtype: {original_dtype}")
print(f"  Converted dtype: {dgrid.dtype}")
```

### Result Extraction Pattern

```python
# Run periodic segmentation
results = mgr.do_periodic_segmentation_rust(
    dgrid, mat_atoms, spacing=spacing,
    mask_thickness=mask_thickness, h=h, threshold_abs=threshold_abs,
    minimum_window_separation=minimum_window_separation,
    wall_thickness=wall_thickness, extra_padding=extra_padding,
    return_padded_regions=True, compute_connections=True
)

# Extract results with clear variable names
maxima_pbc = results['maxima_pbc']           # One per PBC group
regions_pbc = results['regions_pbc']         # Segmentation WITH PBC merged
regions = results['regions']                 # Segmentation WITHOUT PBC (split)
connections = results['connections']         # Windows between pores
all_maxima = results['all_maxima']           # All periodic images
unique_regions = results['unique_regions']   # Regions influencing unit cell
pbc_groups = results['pbc_groups']           # Lists of periodic image groups

# Clever indexing: Maxima influencing the unit cell
maxima = all_maxima[:regions.max()]
```

### Connection Trimming

```python
# Trim to connections that influence the unit cell
connections_trimmed = np.array([
    c for c in connections
    if np.any(np.isin(c[0], unique_regions))
])

print(f"Original connections: {len(connections)}")
print(f"Trimmed connections: {len(connections_trimmed)}")
```

**Why:**
- Removes artifacts from padded regions
- Ensures analysis focuses on relevant pores
- Prevents fictitious results

### Visualization Pattern

```python
import plotly.graph_objects as go
import plotly.express as px

# Set up colorscale
cscale = px.colors.qualitative.Pastel

# Update structure with covalent bonds
mat_updated = mgr.update_with_covalent_bonds(mat_atoms)

# Create base traces
cell_trace = mgr.plot_cell(mat_atoms)
frame_trace = mgr.plot_framework_as_mesh(mat_updated, stylize=True)

# Plot merged regions
surface_trace = mgr.plot_surface_dgrid(
    dgrid,
    isovalue=1.5,
    region_labels=regions_pbc,
    ase_atoms=mat_atoms,
    cscale=cscale,
    name='Pores (PBC merged)',
    stylize=True,
    opacity=0.8
)

# Plot connections (use this, not plot_rag)
connections_traces = mgr.plot_connections_as_it_is(
    connections_trimmed,
    all_maxima,
    dgrid.shape,
    mat_atoms,
    cscale=cscale,
    plot_disconnected=False,
    show_spheres=False
)

# Combine into figure
fig = go.Figure(
    data=[cell_trace, surface_trace, *connections_traces, *frame_trace],
    layout=mgr.dark_layout  # or mgr.white_layout for documentation
)

fig.update_layout(
    scene_aspectmode='data',
    title=dict(
        text="IRMOF-1: Pore Network with Connections",
        x=0.5,
        xanchor='center',
        font=dict(size=16, color='white')
    ),
    scene=dict(
        xaxis=dict(title='X (Å)'),
        yaxis=dict(title='Y (Å)'),
        zaxis=dict(title='Z (Å)')
    )
)

# Save output
fig.write_html('output.html')
fig.write_image('output.png', width=1200, height=900, scale=2)
```

**Key Points:**
- Use `plot_connections_as_it_is()` instead of `plot_rag()` for better results
- Use `dark_layout` for interactive viewing, `white_layout` for documentation
- Always save both HTML (interactive) and PNG (static)

---

## Image Generation

### Why Generate Images?

**Images are CRITICAL for tutorials:**
- Show expected outputs
- Break up text
- Provide concrete examples
- Build confidence ("this is what it should look like")

### Image Generation Script Pattern

Create a separate script (`generate_tutorial_images.py`) that:

1. **Runs the full workflow**
2. **Saves images with white backgrounds** (for documentation)
3. **Uses descriptive filenames**
4. **Saves to a dedicated directory** (`tutorial_images/`)

**Example structure:**

```python
#!/usr/bin/env python3
"""
Generate images for the [Tutorial Name] tutorial.

This script creates PNG visualizations using the workflow.
Images are saved with white backgrounds for documentation use.
"""

import os
import sys
import numpy as np
from pathlib import Path

# Setup paths
script_path = Path(__file__).resolve()
if script_path.parent.name == 'docs':
    repo_root = script_path.parent.parent
else:
    repo_root = script_path.parent

sys.path.insert(0, str(repo_root))

import mofography as mgr
from ase.io import read
import plotly.graph_objects as go
import plotly.express as px

# =============================================================================
# Load Structure
# =============================================================================

print("Loading structure...")
path_to_cif = repo_root / 'examples' / 'TestCIFs' / 'P1-IRMOF-1.cif'
mat_atoms = read(str(path_to_cif))

# =============================================================================
# Run Workflow
# =============================================================================

print("\nRunning workflow...")

# [Workflow code here - same as tutorial]

# =============================================================================
# Create output directory
# =============================================================================

output_dir = Path(__file__).parent / 'tutorial_images'
output_dir.mkdir(exist_ok=True)
print(f"\nSaving images to: {output_dir}")

# =============================================================================
# Generate Image 1: [Description]
# =============================================================================

print("\nGenerating image: image1.png")

# [Create visualization]

fig1 = go.Figure(
    data=[...],
    layout=mgr.white_layout  # Use white layout for documentation
)

fig1.update_layout(
    scene_aspectmode='data',
    title=dict(
        text="[Descriptive Title]",
        x=0.5,
        xanchor='center',
        font=dict(size=14, color='black')
    )
)

fig1.write_image(output_dir / 'image1.png', width=800, height=600, scale=2)
print("  Saved: image1.png")

# [Repeat for other images]

print("\n✅ All images generated successfully!")
print(f"Images saved to: {output_dir}")
```

### Image Best Practices

1. **Use white backgrounds** (`mgr.white_layout`) for documentation
2. **High resolution**: `width=800, height=600, scale=2` or higher
3. **Descriptive filenames**: `split_regions.png`, `merged_regions.png`
4. **Descriptive titles**: Each figure should have a clear title
5. **Consistent styling**: Use same colorscale across related images
6. **Show progression**: Before/after images are powerful

### Including Images in Markdown

```markdown
![Descriptive alt text](./tutorial_images/image_name.png)

*Caption explaining what the image shows and why it matters*
```

**Important:**
- Use relative paths from the markdown file location
- Alt text should be descriptive
- Include captions for context

---

## Common Pitfalls

### 1. Images Not Showing

**Problem**: Images don't render in markdown

**Solutions:**
- Check relative paths are correct from markdown file location
- Ensure images actually exist in the specified directory
- Use `./` for current directory references
- Verify filename spelling matches exactly (case-sensitive!)

### 2. Wrong Dtype

**Problem**: Rust backend expects float32, gets float64

**Solution**:
```python
dgrid = dgrid.astype(np.float32)  # CRITICAL step
```

### 3. Path Issues

**Problem**: CIF files not found when running from different directories

**Solution**:
- Use `pathlib.Path` instead of `os.path`
- Always use `.expanduser()` for tilde paths
- Support multiple possible locations
- Allow environment variable override

### 4. Missing Dependencies

**Problem**: Code fails due to missing packages

**Solution**: Include prerequisites section with:
- Required packages
- Installation instructions
- Version information if critical

### 5. Too Much Information

**Problem**: Tutorial tries to cover everything

**Solution**:
- Use "Tutorial Scope" section to set boundaries
- Link to more detailed resources
- Focus on "quick start" not "exhaustive reference"

### 6. Dry, Boring Content

**Problem**: Tutorial is technically correct but uninspiring

**Solution**:
- Add "Why this matters" sections
- Use conversational tone
- Include emojis and visual breaks
- Show, don't just tell (with images!)
- Use concrete examples

### 7. Code That Doesn't Run

**Problem**: Copy-pasted code fails

**Solution**:
- Test ALL code before including
- Include all imports
- Define all variables
- Use complete, runnable examples
- Test in clean environment

### 8. Wrong Visualization Functions

**Problem**: Using `plot_rag()` when `plot_connections_as_it_is()` is better

**Solution**:
- Test visualization approaches
- Use the function that gives clearest results
- Document why you chose a particular approach

---

## Validation Checklist

Before publishing a tutorial, verify:

### Content
- [ ] Title is catchy and descriptive
- [ ] Introduction explains WHY before HOW
- [ ] Scope is clearly defined
- [ ] All code is runnable as-is
- [ ] All parameters are explained
- [ ] Troubleshooting section included
- [ ] Next steps section provides resources

### Code
- [ ] All imports included
- [ ] All variables defined before use
- [ ] Code tested in clean environment
- [ ] Critical type conversions documented
- [ ] Error handling where appropriate
- [ ] Output expectations shown

### Images
- [ ] All images exist and are accessible
- [ ] Image paths are correct relative to markdown file
- [ ] Images have descriptive alt text
- [ ] Images have captions explaining their relevance
- [ ] Image resolution is sufficient (scale=2 recommended)

### Style
- [ ] Consistent emoji usage
- [ ] Code blocks have language specified
- [ ] Inline code uses backticks
- [ ] Tables for structured information
- [ ] Conversational but professional tone
- [ ] No walls of text (break up with headers/visuals)

### Links
- [ ] All internal links work
- [ ] All external links are accessible
- [ ] References to other tutorials are accurate

### Accessibility
- [ ] Images have alt text
- [ ] Color contrasts are sufficient
- [ ] Code is readable (appropriate line length)
- [ ] Explanations assume no prior knowledge

---

## Summary: The Quick Start Template

When creating a new tutorial, start with this structure:

```markdown
# 🌟 [Catchy Title]

> **[One-sentence summary]**

---

## 📚 Citation

[Paper info if applicable]

---

## 🎯 What is [Topic]?

[Core concept explanation]

### Why is this Powerful? 🚀

[Benefits/bullet points]

---

## ⚠️ Tutorial Scope

**This is a [TYPE] tutorial**

✅ Covers:
- [Topics]

❌ Does NOT cover:
- [Topics]

---

## 🔧 Prerequisites

[Requirements]

---

## 📖 Step-by-Step

### Step 1: [Title] 🔧

[What and why]

```python
[code]
```

[What happened / what is this]

---

## 🚀 Next Steps

[Links to related content]

---

## 📝 Summary

[Table format summary]

**Key Takeaway:** [Main point]

---

## 🎁 Bonus: [Title]

[Extra helpful content]
```

---

## Final Advice

1. **Start with WHY** - Explain the problem before the solution
2. **Show, don't tell** - Use images and examples
3. **Keep it focused** - Define scope and stick to it
4. **Test everything** - Code must run, images must load
5. **Be conversational** - Write like you're explaining to a colleague
6. **Use markdown** - Better than Jupyter for tutorials
7. **Generate images** - Make them white-background, high-res
8. **Learn from feedback** - Iterate and improve

---

**Created**: December 2024
**Based on**: Periodic Segmentation Quick Start Tutorial
**Purpose**: Guide for future tutorial creation in MOFography
