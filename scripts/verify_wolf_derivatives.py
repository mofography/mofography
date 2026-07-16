#!/usr/bin/env python3
"""Verify Wolf DSF radial derivative coefficients against SymPy."""

import sympy as sp

r, a, Rc, qq = sp.symbols("r a Rc qq", positive=True, real=True)
er = sp.erfc(a * r)
eR = sp.erfc(a * Rc)
S = eR / Rc
M = S / Rc + 2 * a * sp.exp(-a**2 * Rc**2) / (sp.sqrt(sp.pi) * Rc)
V = qq * (er / r - S + M * (r - Rc))
f1 = sp.simplify(sp.diff(V, r) / r)
f2 = sp.simplify(sp.diff(f1, r) / r)
f3 = sp.simplify(sp.diff(f2, r) / r)

print("f1 =", f1)
print("f2 =", f2)
print("f3 =", f3)

# Spot-check at numeric point
subs = {r: sp.Float(3.5), a: sp.Float(0.32), Rc: sp.Float(12.5), qq: sp.Float(0.42)}
for name, expr in [("f1", f1), ("f2", f2), ("f3", f3)]:
    print(name, float(expr.subs(subs)))
