#!/usr/bin/env ipython

import stackprinter
stackprinter.set_excepthook(style='darkbg2')
import mofography as mgr
from ase.io import read
from pathlib import Path
path_to_cif = Path('~/Datasets/IZA_cifs/ABW-0.cif').expanduser()
atoms = read(path_to_cif)
dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=0.5)
data= mgr.do_periodic_segmentation_rust(dgrid, spacing=0.5, mat_atoms =atoms, h=0.2, mask_thickness=0.5, threshold_abs = 1.0, extra_padding=5, )
sa = mgr.compute_surface_area_from_results(data, atoms, probe_radius = 1.2)
print(sa)
