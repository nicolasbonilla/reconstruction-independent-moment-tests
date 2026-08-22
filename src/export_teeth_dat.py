# -*- coding: utf-8 -*-
"""Write fig_teeth.dat from the sector-Lanczos teeth cache (spectral_lanczos.run_teeth):
coverage (%), independent-estimator first-moment discrepancy (%), and the circular control (0)."""
import os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
OUT = os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs'))
L = int(os.environ.get('TEETH_L', '12'))


def main():
    dat = np.load(os.path.join(CACHE, f'teeth_L{L}.npz'))
    dd, cov, rind, rcirc = dat['d'], dat['cov'], dat['rind'], dat['rcirc']
    m1_op = float(dat['m1_op']); nsup = int(dat['nsup'])
    with open(os.path.join(OUT, 'teeth.dat'), 'w') as f:
        f.write("d cov rind rcirc\n")
        for dv, c, ri, rc in zip(dd, cov, rind, rcirc):
            f.write(f"{int(dv)} {100*c:.4f} {100*ri:.4f} {100*rc:.6f}\n")
    # report the first threshold (5%) crossing subspace size d for the caption/annotation
    cross = None
    for k in range(1, len(dd)):
        if 100 * rind[k - 1] < 5.0 <= 100 * rind[k]:
            import numpy as _np
            lo, hi = _np.log10(dd[k - 1]), _np.log10(dd[k]); r0, r1 = 100 * rind[k - 1], 100 * rind[k]
            cross = 10 ** (lo + (hi - lo) * (5.0 - r0) / (r1 - r0))
            break
    print(f"L={L}  m1_op={m1_op:.4f}  n_support={nsup}  worst r_indep={100*rind[-1]:.1f}%")
    print(f"first 5% crossing at d ~ {cross:.0f}  (cov ~ {100*cross/nsup:.1f}%)" if cross else "no clean 5% crossing")
    print("wrote teeth.dat ->", OUT)


if __name__ == '__main__':
    main()
