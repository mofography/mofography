# Van der Waals Radii Comparison Report

This report compares the Van der Waals (vdW) radii defined in Zeo++ with values from a comparison table provided in `Docs/vdw_radii_comparison.md`.

## 1. Location of Zeo++ Radii
Let's start by pointing out the two sets of radii that ships with zeo++0.3. The default radii and teh `UFF.rad` file that is included. 

*   **Zeo++ Default Radii:** These values are hardcoded in the `initializeRadTable()` function within the `zeo++-0.3/networkinfo.cc` file.
*   **Zeo++ UFF Radii:** These values are read from the `zeo++-0.3/UFF.rad` file.

## 2. Zeo++ Radii Dictionaries (in Angstroms)

### Zeo++ Default Radii
```python
zeo_default_radii = {
    "H": 1.09, "D": 1.09, "He": 1.4, "Li": 1.82, "Be": 2.0, "B": 2.0, "C": 1.7, "N": 1.55, "O": 1.52, "F": 1.47,
    "Ne": 1.54, "Na": 2.27, "Mg": 1.73, "Al": 2.0, "Si": 2.1, "P": 1.8, "S": 1.8, "Cl": 1.75, "Ar": 1.88, "K": 2.75,
    "Ca": 2.0, "Sc": 2.0, "Ti": 2.0, "V": 2.0, "Cr": 2.0, "Mn": 2.0, "Fe": 2.0, "Co": 2.0, "Ni": 1.63, "Cu": 1.4,
    "Zn": 1.39, "Ga": 1.87, "Ge": 2.0, "As": 1.85, "Se": 1.9, "Br": 1.85, "Kr": 2.02, "Rb": 2.0, "Sr": 2.0, "Y": 2.0,
    "Zr": 2.0, "Nb": 2.0, "Mo": 2.0, "Tc": 2.0, "Ru": 2.0, "Rh": 2.0, "Pd": 1.63, "Ag": 1.72, "Cd": 1.58, "In": 1.93,
    "Sn": 2.17, "Sb": 2.0, "Te": 2.06, "I": 1.98, "Xe": 2.16, "Cs": 2.0, "Ba": 2.0, "La": 2.0, "Ce": 2.0, "Pr": 2.0,
    "Nd": 2.0, "Pm": 2.0, "Sm": 2.0, "Eu": 2.0, "Gd": 2.0, "Tb": 2.0, "Dy": 2.0, "Ho": 2.0, "Er": 2.0, "Tm": 2.0,
    "Yb": 2.0, "Lu": 2.0, "Hf": 2.0, "Ta": 2.0, "W": 2.0, "Re": 2.0, "Os": 2.0, "Ir": 2.0, "Pt": 1.72, "Au": 1.66,
    "Hg": 1.55, "Tl": 1.96, "Pb": 2.02, "Bi": 2.0, "Po": 2.0, "At": 2.0, "Rn": 2.0, "Fr": 2.0, "Ra": 2.0, "Ac": 2.0,
    "Th": 2.0, "Pa": 2.0, "U": 1.86, "Np": 2.0, "Pu": 2.0, "Am": 2.0, "Cm": 2.0, "Bk": 2.0, "Cf": 2.0, "Es": 2.0,
    "Fm": 2.0, "Md": 2.0, "No": 2.0, "Lr": 2.0, "Rf": 2.0, "Db": 2.0, "Sg": 2.0, "Bh": 2.0, "Hs": 2.0, "Mt": 2.0,
    "Ds": 2.0
}
```

### Zeo++ UFF Radii
```python
zeo_uff_radii = {
    "Ac": 1.55, "Ag": 1.4, "Al": 2.005, "Am": 1.505, "Ar": 1.725, "As": 1.885, "At": 2.115, "Au": 1.465, "B": 1.82,
    "Ba": 1.65, "Be": 1.225, "Bi": 1.945, "Bk": 1.485, "Br": 1.865, "C": 1.715, "Ca": 1.515, "Cd": 1.27, "Ce": 1.585,
    "Cf": 1.475, "Cl": 1.76, "Cm": 1.48, "Co": 1.28, "Cr": 1.345, "Cu": 1.555, "Cs": 2.01, "Dy": 1.525, "Eu": 1.555,
    "Er": 1.51, "Es": 1.47, "F": 1.5, "Fe": 1.295, "Fm": 1.465, "Fr": 2.185, "Ga": 1.95, "Ge": 1.905, "Gd": 1.5,
    "H": 1.285, "He": 1.32, "Hf": 1.4, "Hg": 1.205, "Ho": 1.52, "I": 2.005, "In": 1.99, "Ir": 1.265, "K": 1.7,
    "Kr": 1.845, "La": 1.57, "Li": 1.09, "Lu": 1.62, "Lr": 1.44, "Md": 1.46, "Mg": 1.345, "Mn": 1.32, "Mo": 1.36,
    "N": 1.63, "Na": 1.33, "Ne": 1.33, "Nb": 1.41, "Nd": 1.59, "No": 1.445, "Ni": 1.26, "Np": 1.525, "O": 1.56,
    "Os": 1.39, "P": 1.845, "Pa": 1.525, "Pb": 1.915, "Pd": 1.29, "Pm": 1.58, "Po": 2.1, "Pr": 1.605, "Pt": 1.225,
    "Pu": 1.525, "Ra": 1.64, "Rb": 1.835, "Re": 1.315, "Rh": 1.305, "Rn": 2.125, "Ru": 1.32, "S": 1.795, "Sb": 1.97,
    "Sc": 1.47, "Se": 1.875, "Si": 1.915, "Sm": 1.57, "Sn": 1.955, "Sr": 1.62, "Ta": 1.41, "Tb": 1.535, "Tc": 1.335,
    "Te": 1.99, "Th": 1.515, "Ti": 1.415, "Tl": 1.935, "Tm": 1.505, "U": 1.51, "V": 1.4, "W": 1.365, "Xe": 1.96,
    "Y": 1.49, "Yb": 1.495, "Zn": 1.23, "Zr": 1.39
}
```

## 3. Comparison of Van der Waals Radii

This section compares the Zeo++ radii with values from `compute.py`, Bondi, Hu, Alvarez, UFF Sigma, and UFF vdW Radius (from `Docs/vdw_radii_comparison.md`).

NOTE: The UFF vdw is computed as $\text{UFF}~\sigma \cdot \dfrac{2^{1/6}}{2}$

| Element | Zeo++ Default | Zeo++ UFF | compute.py | Bondi | Hu | Alvarez | UFF Sigma | UFF vdW Radius |
|---|---|---|---|---|---|---|---|---|
| H | 1.09 | 1.285 | 1.10 | 1.2 | 1.08 | 1.2 | 2.57113 | 1.447 |
| He | 1.4 | 1.32 | 1.40 |  |  |  | 2.10430 | 1.185 |
| Li | 1.82 | 1.09 | 1.82 | 1.82 |  |  |  |  |
| Be | 2.0 | 1.225 | 1.53 |  |  |  | 2.44552 | 1.375 |
| B | 2.0 | 1.82 | 1.92 |  |  |  | 3.63754 | 2.046 |
| C | 1.7 | 1.715 | 1.70 | 1.7 | 1.53 | 1.77 | 3.43085 | 1.931 |
| N | 1.55 | 1.63 | 1.55 | 1.55 | 1.51 | 1.66 | 3.26069 | 1.835 |
| O | 1.52 | 1.56 | 1.52 | 1.52 | 1.49 | 1.5 | 3.11815 | 1.755 |
| F | 1.47 | 1.5 | 1.47 | 1.47 | 1.48 | 1.46 | 2.99698 | 1.687 |
| Ne | 1.54 | 1.33 | 1.54 |  |  |  |  |  |
| Na | 2.27 | 1.33 | 2.27 | 2.27 |  |  | 2.65755 | 1.494 |
| Mg | 1.73 | 1.345 | 1.73 | 1.73 |  |  | 2.69141 | 1.513 |
| Al | 2.0 | 2.005 | 1.84 | 2.05 |  |  | 4.00815 | 2.257 |
| Si | 2.1 | 1.915 | 2.10 | 2.1 | 1.93 | 2.19 | 3.82641 | 2.154 |
| P | 1.8 | 1.845 | 1.80 | 1.8 | 1.88 | 1.9 | 3.69456 | 2.078 |
| S | 1.8 | 1.795 | 1.80 | 1.8 | 1.81 | 1.89 | 3.59478 | 2.022 |
| Cl | 1.75 | 1.76 | 1.75 | 1.75 | 1.75 | 1.82 | 3.51638 | 1.978 |
| Ar | 1.88 | 1.725 | 1.88 |  |  |  | 3.44600 | 1.938 |
| K | 2.75 | 1.7 | 2.75 | 2.75 |  |  | 3.39611 | 1.910 |
| Ca | 2.0 | 1.515 | 2.31 | 2.45 |  |  | 3.02816 | 1.703 |
| Sc | 2.0 | 1.47 | 2.15 | 2.25 |  |  | 2.93551 | 1.650 |
| Ti | 2.0 | 1.415 | 2.11 | 2.10 |  |  | 2.82860 | 1.591 |
| V | 2.0 | 1.4 | 2.07 | 2.05 |  |  | 2.80099 | 1.575 |
| Cr | 2.0 | 1.345 | 2.06 | 2.0 |  |  | 2.69319 | 1.512 |
| Mn | 2.0 | 1.32 | 2.05 | 2.0 |  |  | 2.63795 | 1.484 |
| Fe | 2.0 | 1.295 | 2.04 | 2.0 |  |  | 2.59430 | 1.459 |
| Co | 2.0 | 1.28 | 2.00 |  |  |  | 2.55866 | 1.441 |
| Ni | 1.63 | 1.26 | 1.97 | 1.63 |  | 2.4 | 2.52481 | 1.419 |
| Cu | 1.4 | 1.555 | 1.96 | 1.40 |  | 2.38 | 3.11369 | 1.753 |
| Zn | 1.39 | 1.23 | 2.01 | 1.39 |  | 2.39 | 2.46155 | 1.386 |
| Ga | 1.87 | 1.95 | 1.87 | 1.87 |  |  | 3.90481 | 2.195 |
| Ge | 2.0 | 1.905 | 2.11 | 2.05 |  |  | 3.81305 | 2.144 |
| As | 1.85 | 1.885 | 1.85 | 1.85 |  | 1.88 | 3.76850 | 2.122 |
| Se | 1.9 | 1.875 | 1.90 | 1.90 |  |  | 3.74623 | 2.109 |
| Br | 1.85 | 1.865 | 1.85 | 1.85 |  | 1.86 | 3.73197 | 2.101 |
| Kr | 2.02 | 1.845 | 2.02 |  |  |  | 3.68921 | 2.077 |
| Rb | 2.0 | 1.835 | 3.03 | 3.0 |  |  | 3.66516 | 2.064 |
| Sr | 2.0 | 1.62 | 2.49 | 2.6 |  |  | 3.24376 | 1.824 |
| Y | 2.0 | 1.49 | 2.32 | 2.4 |  |  | 2.98006 | 1.676 |
| Zr | 2.0 | 1.39 | 2.23 | 2.3 |  |  | 2.78317 | 1.565 |
| Nb | 2.0 | 1.41 | 2.18 | 2.15 |  |  | 2.81969 | 1.585 |
| Mo | 2.0 | 1.36 | 2.17 | 2.1 |  |  | 2.71902 | 1.530 |
| Tc | 2.0 | 1.335 | 2.16 | 2.1 |  |  | 2.67091 | 1.502 |
| Ru | 2.0 | 1.32 | 2.13 | 2.05 |  |  | 2.63973 | 1.484 |
| Rh | 2.0 | 1.305 | 2.10 | 2.0 |  |  | 2.60944 | 1.469 |
| Pd | 1.63 | 1.29 | 2.10 | 1.63 |  | 2.15 | 2.58272 | 1.452 |
| Ag | 1.72 | 1.4 | 2.11 | 1.72 |  | 2.53 | 2.80455 | 1.577 |
| Cd | 1.58 | 1.27 | 2.18 | 1.62 |  | 2.43 | 2.53728 | 1.427 |
| In | 1.93 | 1.99 | 1.93 |  |  | 2.43 | 3.97608 | 2.240 |
| Sn | 2.17 | 1.955 | 2.17 |  |  | 2.42 | 3.91283 | 2.195 |
| Sb | 2.0 | 1.97 | 2.06 |  |  |  | 3.93777 | 2.214 |
| Te | 2.06 | 1.99 | 2.06 | 2.06 |  | 1.99 | 3.98232 | 2.244 |
| I | 1.98 | 2.005 | 1.98 | 1.98 |  | 2.04 | 4.00904 | 2.257 |
| Xe | 2.16 | 1.96 | 2.16 |  |  |  | 3.92352 | 2.202 |
| Cs | 2.0 | 2.01 | 3.43 | 3.15 |  |  | 4.02419 | 2.263 |
| Ba | 2.0 | 1.65 | 2.68 | 2.7 |  |  | 3.29900 | 1.855 |
| La | 2.0 | 1.57 | 2.43 |  |  |  | 3.13775 | 1.765 |
| Ce | 2.0 | 1.585 | 2.42 |  |  |  | 3.16804 | 1.781 |
| Pr | 2.0 | 1.605 | 2.40 |  |  |  | 3.21258 | 1.806 |
| Nd | 2.0 | 1.59 | 2.39 |  |  |  | 3.18496 | 1.793 |
| Pm | 2.0 | 1.58 | 2.38 |  |  |  | 3.16002 | 1.779 |
| Sm | 2.0 | 1.57 | 2.36 |  |  |  | 3.13596 | 1.765 |
| Eu | 2.0 | 1.555 | 2.35 |  |  |  | 3.11191 | 1.752 |
| Gd | 2.0 | 1.5 | 2.34 |  |  |  | 3.00055 | 1.687 |
| Tb | 2.0 | 1.535 | 2.33 |  |  |  | 3.07449 | 1.732 |
| Dy | 2.0 | 1.525 | 2.31 |  |  |  | 3.05400 | 1.721 |
| Ho | 2.0 | 1.52 | 2.30 |  |  |  | 3.03707 | 1.711 |
| Er | 2.0 | 1.51 | 2.29 |  |  |  | 3.02104 | 1.703 |
| Tm | 2.0 | 1.505 | 2.27 |  |  |  | 3.00589 | 1.695 |
| Yb | 2.0 | 1.495 | 2.26 |  |  |  | 2.98897 | 1.685 |
| Lu | 2.0 | 1.62 | 2.24 |  |  |  | 3.24287 | 1.824 |
| Hf | 2.0 | 1.4 | 2.23 |  |  |  | 2.79831 | 1.574 |
| Ta | 2.0 | 1.41 | 2.22 |  |  |  | 2.82415 | 1.588 |
| W | 2.0 | 1.365 | 2.18 |  |  |  | 2.73417 | 1.538 |
| Re | 2.0 | 1.315 | 2.16 |  |  |  | 2.63171 | 1.481 |
| Os | 2.0 | 1.39 | 2.16 |  |  |  | 2.77960 | 1.563 |
| Ir | 2.0 | 1.265 | 2.13 |  |  |  | 2.53015 | 1.424 |
| Pt | 1.72 | 1.225 | 2.13 |  |  |  | 2.45354 | 1.382 |
| Au | 1.66 | 1.465 | 2.14 |  |  |  | 2.93373 | 1.650 |
| Hg | 1.55 | 1.205 | 2.23 |  |  |  | 2.40988 | 1.357 |
| Tl | 1.96 | 1.935 | 1.96 |  |  |  | 3.87274 | 2.181 |
| Pb | 2.02 | 1.915 | 2.02 |  |  |  | 3.82819 | 2.154 |
| Bi | 2.0 | 1.945 | 2.07 |  |  |  | 3.89323 | 2.188 |
| Po | 2.0 | 2.1 | 1.97 |  |  |  | 4.19524 | 2.362 |
| At | 2.0 | 2.115 | 2.02 |  |  |  | 4.23177 | 2.382 |
| Rn | 2.0 | 2.125 | 2.20 |  |  |  | 4.24513 | 2.389 |
| Fr | 2.0 | 2.185 | 3.48 |  |  |  | 4.36540 | 2.459 |
| Ra | 2.0 | 1.64 | 2.83 |  |  |  |  |  |
| Ac | 2.0 | 1.55 | 2.47 |  |  |  | 3.09855 | 1.742 |
| Th | 2.0 | 1.515 | 2.45 |  |  |  | 3.02549 | 1.703 |
| Pa | 2.0 | 1.525 | 2.43 |  |  |  | 3.05044 | 1.717 |
| U | 1.86 | 1.51 | 2.41 |  |  |  | 3.02460 | 1.703 |
| Np | 2.0 | 1.525 | 2.39 |  |  |  | 3.05044 | 1.717 |
| Pu | 2.0 | 1.525 | 2.43 |  |  |  | 3.05044 | 1.717 |
| Am | 2.0 | 1.505 | 2.44 |  |  |  | 3.01213 | 1.695 |
| Cm | 2.0 | 1.48 | 2.45 |  |  |  | 2.96313 | 1.667 |

## 4. Notes on Significant Discrepancies

This section highlights elements where the Van der Waals radii in Zeo++ differ significantly (by more than 0.1 Å) from the reference values in the comparison table.

* **Ac**: Zeo++ Default (2.0) vs compute.py (2.47)
* **Ag**: Zeo++ Default (1.72) vs Alvarez (2.53)
* **Ag**: Zeo++ Default (1.72) vs compute.py (2.11)
* **Al**: Zeo++ Default (2.0) vs compute.py (1.84)
* **Am**: Zeo++ Default (2.0) vs compute.py (2.44)
* **Ar**: Zeo++ UFF (1.725) vs UFF vdW Radius (1.938)
* **Au**: Zeo++ Default (1.66) vs compute.py (2.14)
* **Ba**: Zeo++ Default (2.0) vs Bondi (2.7)
* **Ba**: Zeo++ Default (2.0) vs compute.py (2.68)
* **Be**: Zeo++ Default (2.0) vs compute.py (1.53)
* **Ca**: Zeo++ Default (2.0) vs Bondi (2.45)
* **Ca**: Zeo++ Default (2.0) vs compute.py (2.31)
* **Cd**: Zeo++ Default (1.58) vs Alvarez (2.43)
* **Cd**: Zeo++ Default (1.58) vs compute.py (2.18)
* **Ce**: Zeo++ Default (2.0) vs compute.py (2.42)
* **Cm**: Zeo++ Default (2.0) vs compute.py (2.45)
* **Co**: Zeo++ UFF (1.28) vs UFF vdW Radius (1.441)
* **Cr**: Zeo++ UFF (1.345) vs UFF vdW Radius (1.512)
* **Cs**: Zeo++ Default (2.0) vs Bondi (3.15)
* **Cs**: Zeo++ Default (2.0) vs compute.py (3.43)
* **Cu**: Zeo++ Default (1.4) vs Alvarez (2.38)
* **Cu**: Zeo++ Default (1.4) vs compute.py (1.96)
* **Dy**: Zeo++ Default (2.0) vs compute.py (2.31)
* **Er**: Zeo++ Default (2.0) vs compute.py (2.29)
* **Eu**: Zeo++ Default (2.0) vs compute.py (2.35)
* **Fe**: Zeo++ UFF (1.295) vs UFF vdW Radius (1.459)
* **Fr**: Zeo++ Default (2.0) vs compute.py (3.48)
* **Gd**: Zeo++ Default (2.0) vs compute.py (2.34)
* **Ge**: Zeo++ Default (2.0) vs compute.py (2.11)
* **H**: Zeo++ UFF (1.285) vs UFF vdW Radius (1.447)
* **Hf**: Zeo++ Default (2.0) vs compute.py (2.23)
* **Hg**: Zeo++ Default (1.55) vs compute.py (2.23)
* **Ho**: Zeo++ Default (2.0) vs compute.py (2.30)
* **In**: Zeo++ Default (1.93) vs Alvarez (2.43)
* **Ir**: Zeo++ UFF (1.265) vs UFF vdW Radius (1.424)
* **K**: Zeo++ UFF (1.7) vs UFF vdW Radius (1.910)
* **La**: Zeo++ Default (2.0) vs compute.py (2.43)
* **Li**: Zeo++ UFF (1.09) vs UFF vdW Radius ()
* **Lu**: Zeo++ Default (2.0) vs compute.py (2.24)
* **Mg**: Zeo++ UFF (1.345) vs UFF vdW Radius (1.513)
* **Mn**: Zeo++ UFF (1.32) vs UFF vdW Radius (1.484)
* **Mo**: Zeo++ Default (2.0) vs compute.py (2.17)
* **Na**: Zeo++ UFF (1.33) vs UFF vdW Radius (1.494)
* **Nb**: Zeo++ Default (2.0) vs Bondi (2.15)
* **Nb**: Zeo++ Default (2.0) vs compute.py (2.18)
* **Nd**: Zeo++ Default (2.0) vs compute.py (2.39)
* **Ni**: Zeo++ Default (1.63) vs Alvarez (2.4)
* **Ni**: Zeo++ Default (1.63) vs compute.py (1.97)
* **Ni**: Zeo++ UFF (1.26) vs UFF vdW Radius (1.419)
* **Np**: Zeo++ Default (2.0) vs compute.py (2.39)
* **Pa**: Zeo++ Default (2.0) vs compute.py (2.43)
* **Pd**: Zeo++ Default (1.63) vs Alvarez (2.15)
* **Pd**: Zeo++ Default (1.63) vs compute.py (2.10)
* **Pm**: Zeo++ Default (2.0) vs compute.py (2.38)
* **Pr**: Zeo++ Default (2.0) vs compute.py (2.40)
* **Pt**: Zeo++ Default (1.72) vs compute.py (2.13)
* **Pu**: Zeo++ Default (2.0) vs compute.py (2.43)
* **Ra**: Zeo++ Default (2.0) vs compute.py (2.83)
* **Rb**: Zeo++ Default (2.0) vs Bondi (3.0)
* **Rb**: Zeo++ Default (2.0) vs compute.py (3.03)
* **Re**: Zeo++ Default (2.0) vs compute.py (2.16)
* **Rh**: Zeo++ Default (2.0) vs compute.py (2.10)
* **Ru**: Zeo++ Default (2.0) vs compute.py (2.13)
* **Sc**: Zeo++ Default (2.0) vs Bondi (2.25)
* **Sc**: Zeo++ Default (2.0) vs compute.py (2.15)
* **Sm**: Zeo++ Default (2.0) vs compute.py (2.36)
* **Sn**: Zeo++ Default (2.17) vs Alvarez (2.42)
* **Sr**: Zeo++ Default (2.0) vs Bondi (2.6)
* **Sr**: Zeo++ Default (2.0) vs compute.py (2.49)
* **Ta**: Zeo++ Default (2.0) vs compute.py (2.22)
* **Tb**: Zeo++ Default (2.0) vs compute.py (2.33)
* **Tc**: Zeo++ Default (2.0) vs Bondi (2.1)
* **Tc**: Zeo++ Default (2.0) vs compute.py (2.16)
* **Th**: Zeo++ Default (2.0) vs compute.py (2.45)
* **Ti**: Zeo++ Default (2.0) vs Bondi (2.10)
* **Ti**: Zeo++ Default (2.0) vs compute.py (2.11)
* **Tm**: Zeo++ Default (2.0) vs compute.py (2.27)
* **U**: Zeo++ Default (1.86) vs compute.py (2.41)
* **W**: Zeo++ Default (2.0) vs compute.py (2.18)
* **Y**: Zeo++ Default (2.0) vs Bondi (2.4)
* **Y**: Zeo++ Default (2.0) vs compute.py (2.32)
* **Yb**: Zeo++ Default (2.0) vs compute.py (2.26)
* **Zn**: Zeo++ Default (1.39) vs Alvarez (2.39)
* **Zn**: Zeo++ Default (1.39) vs compute.py (2.01)
* **Zr**: Zeo++ Default (2.0) vs Bondi (2.3)
* **Zr**: Zeo++ Default (2.0) vs compute.py (2.23)
