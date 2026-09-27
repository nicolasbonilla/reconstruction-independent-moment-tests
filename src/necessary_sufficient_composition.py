r"""
R5 (2026-09-27 revision): Hankel lower bound and Hausdorff (interval-localizing) upper bound on m_2 for the
L=6, U/t=8 doped-current Lehmann measure -- TWO NECESSARY conditions on one moment sequence.

The 2026-08 version of this script (kept verbatim in src/_superseded/necessary_sufficient_composition_2026-08.py;
its output data/necessary_sufficient_composition.json is kept unchanged as the record) called the upper bound the
"sufficient (Wang/Mortimer SDP)" layer and took the support [a,b] from ALL positive-frequency eigenvalues of the
full Fock space, including zero-weight states up to 52.21 t, giving m2 <= 573.68. This revision:
  * restates the role: for a positive measure on [a,b] with moments m0,m1,m2,
        m1^2/m0 <= m2 <= (a+b) m1 - a b m0          (Hankel PSD; E[(w-a)(b-w)] >= 0)
    both are NECESSARY conditions (a moment vector violating either is not the moment vector of any positive
    measure on [a,b]); neither certifies a reconstruction, and this is not the Wang et al. / Mortimer et al.
    semidefinite programme, whose composition with the screen is not implemented;
  * evaluates the upper bound for four support intervals, from a-priori to oracle:
        (i)  [0, E_max - E0] full Fock range (valid a priori for any T=0 spectrum of this H);
        (ii) [w_min>0, E_max - E0] of the full Fock space (the 2026-08 choice; reproduces 573.68);
        (iii)[0, max excitation of the (N=4, N_up=2) symmetry sector that J|0> lives in] (a priori from symmetry);
        (iv) the WEIGHTED support {w_n > 1e-12} = [2.49, 13.92] t (oracle: needs the true spectrum) -> 147.6;
  * records the weight gap that makes the 1e-12 filter unambiguous.
Writes key 'R5_hausdorff_weighted_support' of data/2026-09-27_theory_numerics.json (never the old JSON).
Deterministic (dense eigensolvers; current_measure() uses ARPACK only to pick the ground state, and the
aggregated measure is M-independent within the spin-triplet ground manifold). Run: cd src && python necessary_sufficient_composition.py
"""
import os, sys, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from markov_krein_window import current_measure   # exact L=6, U/t=8 doped current Lehmann measure (m0=1)
from small_checks import record_theory_numerics

WTHR = 1e-12


def hausdorff_hi(a, b, m0, m1):
    return (a + b) * m1 - a * b * m0


def sector_range():
    """max excitation energy (and lowest positive one) of the (N=4, N_up=2) sector, weight-independent."""
    import hubbard_ed as H
    L, U = 6, 8.0
    c = H.build_operators(L); cd = [ci.getH() for ci in c]
    Ham = H.hubbard_hamiltonian(L, t=1.0, U=U, pbc=True, c=c).tocsr()
    Nd = np.real(sum((cd[q] @ c[q]).diagonal() for q in range(2 * L)))
    Nu = np.real(sum((cd[q] @ c[q]).diagonal() for q in range(0, 2 * L, 2)))
    idx = np.where((np.abs(Nd - 4) < 1e-9) & (np.abs(Nu - 2) < 1e-9))[0]
    E = np.linalg.eigvalsh(np.real(Ham[idx][:, idx].toarray()))
    om = E - E[0]
    return float(om[om > 1e-6].min()), float(om.max()), int(len(idx))


def main():
    t0 = time.time()
    om, w = current_measure()                          # unique poles (omega>1e-6), weights normalized (m0=1)
    m0 = float(np.sum(w)); m1 = float(np.sum(w * om)); m2 = float(np.sum(w * om ** 2))
    hankel_lo = m1 ** 2 / m0
    keep = w > WTHR
    a_w, b_w = float(om[keep].min()), float(om[keep].max())
    a_f, b_f = float(om.min()), float(om.max())
    s_lo, s_hi, s_dim = sector_range()
    supports = {
        'i_full_fock_0_to_Emax': (0.0, b_f),
        'ii_full_fock_min_positive_to_Emax_(2026-08 choice)': (a_f, b_f),
        'iii_symmetry_sector_0_to_max': (0.0, s_hi),
        'iv_weighted_support_oracle': (a_w, b_w),
    }
    bounds = {k: dict(a=a, b=b, m2_upper=hausdorff_hi(a, b, m0, m1),
                      m2_upper_over_true=hausdorff_hi(a, b, m0, m1) / m2,
                      width_of_interval_on_m2=hausdorff_hi(a, b, m0, m1) - hankel_lo)
              for k, (a, b) in supports.items()}
    old = json.load(open(os.path.join(HERE, '..', 'data', 'necessary_sufficient_composition.json')))
    hi_w = bounds['iv_weighted_support_oracle']['m2_upper']
    width = hi_w - hankel_lo
    m2_under = hankel_lo - 0.05 * width
    m2_over = hi_w + 0.05 * width
    rt = time.time() - t0
    R5 = dict(
        measure='L=6, U/t=8 doped current Lehmann measure (markov_krein_window.current_measure), m0 = 1',
        n_atoms_all=int(len(w)), n_atoms_weighted=int(keep.sum()), weight_threshold=WTHR,
        smallest_kept_weight=float(w[keep].min()), largest_dropped_weight=float(w[~keep].max()) if (~keep).any() else 0.0,
        m0=m0, m1=m1, m2=m2, mean=m1 / m0, variance=m2 / m0 - (m1 / m0) ** 2,
        hankel_lower_bound_m2=hankel_lo, true_m2_minus_hankel=m2 - hankel_lo,
        hausdorff_upper_bounds=bounds,
        symmetry_sector=dict(N=4, N_up=2, dim=s_dim, lowest_positive_excitation=s_lo, max_excitation=s_hi),
        reproduces_old_json=dict(old_support=old['support'], old_support_upper_bound_m2=old['support_upper_bound_m2'],
                                 new_ii_value=bounds['ii_full_fock_min_positive_to_Emax_(2026-08 choice)']['m2_upper'],
                                 abs_diff=abs(old['support_upper_bound_m2'] -
                                              bounds['ii_full_fock_min_positive_to_Emax_(2026-08 choice)']['m2_upper'])),
        interval_on_m2_weighted_support=[hankel_lo, hi_w], true_m2_inside=bool(hankel_lo <= m2 <= hi_w),
        demo_corruptions=dict(under_m2=m2_under, under_violates_hankel=bool(m2_under < hankel_lo),
                              over_m2=m2_over, over_violates_hausdorff_weighted=bool(m2_over > hi_w),
                              over_passes_2026_08_bound=bool(m2_over <= old['support_upper_bound_m2'])),
        role='both inequalities are NECESSARY conditions on a moment vector (Hankel PSD; interval/Hausdorff '
             'localizing E[(w-a)(b-w)] >= 0); the combined screen is necessary AND necessary. Not the Wang/Mortimer '
             'sufficient-condition SDP (its composition with the screen is not implemented).',
        caveat='the weighted support [a,b] is ORACLE knowledge (it needs the exact spectrum); an a-priori-valid '
               'interval (i) or (iii) gives a much weaker bound. 147.6 is the tightest Hausdorff bound for this '
               'measure, 573.7/579.7 are the a-priori ones.',
        plan_expectation='[a,b] = [2.49, 13.92] t, m2 <= 147.6 (was 573.7 from zero-weight Fock states up to 52.21 t)')
    print(f"weighted support [{a_w:.4f}, {b_w:.4f}]  ({int(keep.sum())} of {len(w)} atoms; kept min "
          f"{R5['smallest_kept_weight']:.2e}, dropped max {R5['largest_dropped_weight']:.2e})")
    print(f"m0={m0:.6f} m1={m1:.6f} m2={m2:.6f}  Hankel lower {hankel_lo:.4f}")
    for k, v in bounds.items():
        print(f"  {k:55s} [{v['a']:.4f}, {v['b']:.4f}] -> m2 <= {v['m2_upper']:.4f}")
    record_theory_numerics('R5_hausdorff_weighted_support', R5, 'necessary_sufficient_composition.py', rt, seed=None)


if __name__ == '__main__':
    main()
