# -*- coding: utf-8 -*-
"""What the fig_sqw rasters draw below the charge-onset line, and what the frequency windows clip (plan B9).

Reads only src/cache/sqw_L12.npz (the broadened L=12 charge S(q,omega) and spin S^zz(q,omega) on their grids)
and paper/figs/sqw_deltac.dat / sqw_extent.dat, and writes data/2026-09-27_sqw_extent_check.json:

  * the charge onset recomputed with export_sqw_dat.py's rule (lowest omega where max_q S > 5% of the global
    maximum) and checked against sqw_deltac.dat;
  * the largest DISPLAYED colour value below the onset line (display = clip(S/vmax, 0, 1)^0.62 with vmax the
    99.5th percentile, as in export_sqw_dat.py), per frequency band, for the current placement of the rasters
    at their true extents (charge 0-11t, spin 0-2.6t) and for the pre-2026-09-27 placement that stretched them
    to the 10t / 2.4t windows;
  * the fraction of the broadened grid weight that the 10t (charge) and 2.4t (spin) windows clip, overall and
    at the worst q;
  * where the dominant spinon peak of each q is drawn relative to its true frequency under both placements.

The broadening widths eta_c, eta_s are read from the cache. The script separates nothing into poles and
tails: 'below the onset only the broadening tail is visible' rests on the displayed values reported here
being small and rising towards the onset. Committed on 2026-09-27 from the phase-A scratch check (same
numbers). No random numbers. Run from src/:  python check_sqw_extent.py   (under a second)
"""
import os, json, time, platform, sys
import numpy as np

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..'))
GAMMA = 0.62                                   # display gamma, as in export_sqw_dat.py
OUT = os.path.join(ROOT, 'data', '2026-09-27_sqw_extent_check.json')


def main():
    d = np.load(os.path.join(HERE, 'cache', 'sqw_L12.npz'))
    qq, wc, ws, Sc, Ss = d['qq'], d['wc'], d['ws'], d['Sc'], d['Ss']
    prof = Sc.max(axis=0)
    onset = float(wc[np.argmax(prof > 0.05 * prof.max())])
    onset_dat = float(open(os.path.join(ROOT, 'paper', 'figs', 'sqw_deltac.dat')).read().split()[1])
    ext = open(os.path.join(ROOT, 'paper', 'figs', 'sqw_extent.dat')).read().split()
    wc_max, ws_max = float(ext[6]), float(ext[7])
    vmax = float(np.percentile(Sc, 99.5))
    disp = np.clip(Sc / vmax, 0, 1) ** GAMMA
    out = {'_provenance': {'script': 'src/check_sqw_extent.py', 'date': '2026-09-27', 'plan_item': 'B9',
                           'inputs': ['src/cache/sqw_L12.npz', 'paper/figs/sqw_deltac.dat', 'paper/figs/sqw_extent.dat'],
                           'seed': None, 'python': sys.version.split()[0], 'numpy': np.__version__,
                           'platform': platform.platform()},
           'eta_c': float(d['eta_c']), 'eta_s': float(d['eta_s']), 'U': float(d['U']), 'L': int(d['L']),
           'onset_rule': 'lowest omega with max_q S_charge > 5% of its global maximum (export_sqw_dat.py)',
           'onset_recomputed': onset, 'onset_in_sqw_deltac_dat': onset_dat,
           'onset_matches_dat': bool(abs(onset - onset_dat) < 5e-5),
           'extent_true': {'charge_wmax': wc_max, 'spin_wmax': ws_max},
           'display': {'gamma': GAMMA, 'vmax_99p5_percentile': vmax, 'global_max': float(Sc.max())}}
    # rows of the charge raster drawn below the onset line, for both placements
    place = {}
    for name, drawn_top in (('true_extent_0_11t (current)', float(wc[-1])), ('stretched_to_10t (before 2026-09-27)', 10.0)):
        drawn = wc * drawn_top / float(wc[-1])
        below = drawn < onset
        place[name] = {'max_display_value_drawn_below_onset': float(disp[:, below].max()),
                       'true_omega_range_drawn_below_onset': [0.0, float(wc[below].max())]}
    out['charge_below_onset_line'] = place
    bands = []
    pmax = disp.max(axis=0)
    for lo, hi in ((0.0, 4.0), (4.0, 5.0), (5.0, 5.5), (5.5, onset)):
        m = (wc >= lo) & (wc < hi)
        i = int(np.argmax(np.where(m, pmax, -1)))
        bands.append({'omega_range': [lo, hi], 'max_display': float(pmax[m].max()), 'at_omega': float(wc[i]),
                      'at_q': float(qq[int(np.argmax(disp[:, i]))]),
                      'linear_fraction_of_vmax': float(pmax[m].max() ** (1 / GAMMA))})
    out['charge_below_onset_profile_true_extent'] = bands
    # weight clipped by the display windows (broadened grid weight, trapezoid over omega)
    def clipped(S, w, wmax):
        tot = np.trapz(S, w, axis=1); keep = w <= wmax
        cut = tot - np.trapz(S[:, keep], w[keep], axis=1)
        return {'window': wmax, 'fraction_clipped_all_q': float(cut.sum() / tot.sum()),
                'fraction_clipped_worst_q': float((cut / tot).max()), 'worst_q': float(qq[int(np.argmax(cut / tot))])}
    out['window_clipping'] = {'charge': clipped(Sc, wc, 10.0), 'spin': clipped(Ss, ws, 2.4)}
    # spin: dominant peak per q and where it is drawn
    pk = np.array([ws[np.argmax(Ss[i])] for i in range(len(qq))])
    out['spin_peaks'] = {'q': qq.tolist(), 'peak_omega_true': pk.tolist(),
                         'max_offset_drawn_vs_true_stretched_to_2p4t': float(np.max(np.abs(pk - pk * 2.4 / float(ws[-1])))),
                         'max_offset_drawn_vs_true_true_extent': 0.0,
                         'half_pixel_omega_spin': float(ws[-1]) / (len(ws) - 1) / 2,
                         'half_pixel_omega_charge': float(wc[-1]) / (len(wc) - 1) / 2}
    out['_provenance']['runtime_s'] = round(time.time() - T0, 3)
    json.dump(out, open(OUT, 'w'), indent=1)
    cur = place['true_extent_0_11t (current)']
    print(f"onset {onset:.4f} (sqw_deltac.dat {onset_dat:.4f}); eta_c = {out['eta_c']}, eta_s = {out['eta_s']}")
    print(f"current placement: max displayed value below the onset line {cur['max_display_value_drawn_below_onset']:.3f} "
          f"(linear {bands[-1]['linear_fraction_of_vmax']:.3f} of vmax, just below the onset); below 5t at most "
          f"{max(b['max_display'] for b in bands[:2]):.3f}")
    wcx = out['window_clipping']
    print(f"clipped by the windows: charge above 10t {100*wcx['charge']['fraction_clipped_all_q']:.1f}% of the grid weight "
          f"(worst q {100*wcx['charge']['fraction_clipped_worst_q']:.1f}%); spin above 2.4t "
          f"{100*wcx['spin']['fraction_clipped_all_q']:.2f}% (worst q {100*wcx['spin']['fraction_clipped_worst_q']:.2f}%)")
    print('wrote', OUT)


if __name__ == '__main__':
    main()
