r"""
Demonstrated necessary-compose-sufficient composition on ONE shared moment sequence {m_k}.

The paper positions the screen as the NECESSARY-condition layer that composes with the
SUFFICIENT-condition semidefinite certificates of Wang et al. / Mortimer et al., acting
on the same ground-state moment sequence {m_k}. This script turns "would compose" into
"does compose", numerically, WITHOUT overclaiming a joint pipeline: it shows the two
layers act on the identical {m_0,m_1,m_2} and catch COMPLEMENTARY corruption directions.

For a T=0 Lehmann measure supported on [a,b]=[omega_min,omega_max] with moments
m_0,m_1,m_2:
  * NECESSARY (our screen, bare Hankel/Stieltjes positivity):   m_2 >= m_1^2 / m_0
        (the Hankel matrix [[m0,m1],[m1,m2]] >= 0).  Catches UNDER-dispersed corruptions.
  * SUFFICIENT (support-localizing SDP, the Wang/Mortimer feasibility relaxation):
        E[(omega-a)(b-omega)] >= 0  =>  m_2 <= (a+b) m_1 - a b m_0.
        Catches OVER-dispersed corruptions that bare Hankel positivity MISSES.
The certified feasible interval on m_2 given (m_0,m_1) and the support is therefore
  [ m_1^2/m_0 ,  (a+b) m_1 - a b m_0 ],
the exact-moment optimum of the moment SDP (attained by the two-point principal
representation), computed here in closed form for the (m0,m1,m2) case. Both bounds are
functions of the SAME {m_k} -- that is the composition.
"""
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from markov_krein_window import current_measure   # exact L=6, U/t=8 doped current Lehmann measure

om, w = current_measure()                          # unique poles (omega>0), weights (m0=1)
a, b = float(om.min()), float(om.max())            # support endpoints
m0 = float(np.sum(w))
m1 = float(np.sum(w * om))
m2 = float(np.sum(w * om ** 2))

hankel_lo = m1 ** 2 / m0                            # NECESSARY (bare Hankel) lower bound on m2
support_hi = (a + b) * m1 - a * b * m0             # SUFFICIENT (support SDP) upper bound on m2

print("=== necessary-compose-sufficient on the shared {m_0,m_1,m_2} (L=6, U/t=8 current measure) ===")
print(f"support [a,b] = [{a:.4f}, {b:.4f}]")
print(f"true moments : m0={m0:.6f}  m1={m1:.6f}  m2={m2:.6f}")
print(f"NECESSARY (Hankel)  lower bound  m2 >= m1^2/m0          = {hankel_lo:.6f}")
print(f"SUFFICIENT (support) upper bound m2 <= (a+b)m1 - ab m0  = {support_hi:.6f}")
print(f"certified feasible interval for m2 : [{hankel_lo:.6f}, {support_hi:.6f}]")
inside = hankel_lo <= m2 <= support_hi
print(f"true m2 inside certified interval? {inside}  (corroborates, as it must)")

# --- complementary corruptions, both on the SAME {m_k} ---
eps = 0.5 * (support_hi - hankel_lo)
m2_under = hankel_lo - 0.05 * (support_hi - hankel_lo)     # under-dispersed: below Hankel floor
m2_over = support_hi + 0.05 * (support_hi - hankel_lo)     # over-dispersed: above support ceiling


def hankel_psd(m0, m1, m2):
    H = np.array([[m0, m1], [m1, m2]])
    return bool(np.linalg.eigvalsh(H)[0] >= -1e-12)


def support_ok(m0, m1, m2):
    return bool(m2 <= support_hi + 1e-12)


for tag, m2c in [("under-dispersed (m2 below Hankel floor)", m2_under),
                 ("over-dispersed  (m2 above support ceiling)", m2_over)]:
    hk = hankel_psd(m0, m1, m2c); su = support_ok(m0, m1, m2c)
    nec = "PASS" if hk else "REFUTE"
    suf = "PASS" if su else "REFUTE"
    caught_by = "necessary (Hankel)" if not hk else ("sufficient (support SDP)" if not su else "neither")
    print(f"\ncorruption {tag}: m2'={m2c:.4f}")
    print(f"   NECESSARY (our Hankel screen): {nec}   SUFFICIENT (support certificate): {suf}")
    print(f"   -> caught by {caught_by}  (composition on the shared m0,m1,m2)")

print("\nKEY: the two layers act on the IDENTICAL {m_0,m_1,m_2}; bare Hankel positivity catches")
print("under-dispersion, the support-localizing SDP catches over-dispersion. Necessary o sufficient")
print("is thus a demonstrated composition on one moment sequence, not merely a framed one.")

out = os.path.normpath(os.path.join(HERE, "..", "data", "necessary_sufficient_composition.json"))
json.dump({"L": 6, "U": 8.0, "support": [a, b], "m0": m0, "m1": m1, "m2": m2,
           "hankel_lower_bound_m2": hankel_lo, "support_upper_bound_m2": support_hi,
           "certified_interval_m2": [hankel_lo, support_hi], "true_m2_inside": inside,
           "corruption_under_m2": m2_under, "corruption_over_m2": m2_over,
           "under_caught_by_necessary": not hankel_psd(m0, m1, m2_under),
           "over_caught_by_sufficient": not support_ok(m0, m1, m2_over)},
          open(out, "w"), indent=1)
print(f"\nwrote {out}")
