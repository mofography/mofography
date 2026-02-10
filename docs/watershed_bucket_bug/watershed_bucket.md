---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.17.3
  kernelspec:
    display_name: prf
    language: python
    name: python3
---

```python
import pandas as pd
import mofography as mgr
from ase.io import read
from scipy.ndimage import label as scipy_label
import numpy as np
import os
from pathlib import Path
import time
```

```python
import os
os.environ["RAYON_NUM_THREADS"] = "8"
```

```python
#parameter set 2_no smoothing
spacing = 0.1
mask_thickness=0
h= 0.3
threshold_abs= 0.605
preserve_threshold=103.0
smooth_threshold=106.0
max_smooth_radius=2.0
extra_padding= 2
return_dgrid = True
min_distance= 1.0
min_region_size=0.0001
minimum_window_separation=1.5
```

```python
mat_atoms = read("/projects/academic/kaihangs/czheng37/PRF/pore_network/1231/Case4-combined-revised.cif")
radius = pd.Series({"C":172})
radii = (radius[mat_atoms.get_chemical_symbols()].values / 100.0)
dgrid = mgr.dgrid_from_atoms_rust(mat_atoms, spacing=spacing, radii=radii)
results = mgr.do_periodic_segmentation_rust(
    dgrid, 
    mat_atoms,
    mask_thickness=mask_thickness,
    h=h,
    spacing=spacing,
    return_dgrid=False,
    threshold_abs=threshold_abs,
    min_distance=min_distance,
    return_padded_regions=True,
    preserve_threshold=preserve_threshold,
    smooth_threshold=smooth_threshold,
    max_smooth_radius=max_smooth_radius,
    extra_padding=extra_padding,
    compute_connections=True,
    min_region_size=min_region_size,
    minimum_window_separation=minimum_window_separation
)

regions = results["regions"]



regions_old, maxima_old = mgr.regions_from_dgrid(dgrid, mask_thickness,h)
```

```python
np.sum(dgrid<0),np.sum(regions==0),np.sum(regions_old==0)
```
