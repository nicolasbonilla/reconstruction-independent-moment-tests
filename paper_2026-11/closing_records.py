# -*- coding: utf-8 -*-
"""Record checks for the present manuscript (main.tex + supplement.tex of this folder), added on 2026-11-14: numbers
that the text prints in a form the deposit's own check does not anchor, each checked against its JSON key in ../data
with the direction of rounding that its sentence requires:
  'up'     upper bound ("up to", "to a relative X"): the stored value rounded UP at the printed precision;
  'down'   lower end of an enclosure: rounded DOWN;
  'near'   point value: rounded to NEAREST;
  'sci_up' upper bound printed as m x 10^e: mantissa rounded UP at its printed decimals;
  'bool'   a stated fact (value must be True).
Why a separate record set: ../src/verify.py (the deposit's own check, not changed) anchors some of these sentences on
the integer or one-digit forms of the earlier manuscript in ../paper, so in this text those anchors read ABSENT. These
records anchor the new forms. A record whose sentence is missing FAILS (it is this text's own form). The report calls
them closing records; the id prefixes A1, A3, A4 and A5 only group them. Run by verify_prb.py. Reads JSON only."""
import os, re, json, math

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, '..', 'data'))
EPS = 1e-9


def _j(f):
    return json.load(open(os.path.join(DATA, f), encoding='utf-8'))


def _nd(s):
    return len(s.split('.')[1]) if '.' in s else 0


def good(printed, value, rule):
    if rule == 'bool':
        return bool(value)
    if rule == 'sci_up':
        mant, expo = printed
        nd = _nd(mant)
        e = math.floor(math.log10(value))
        m = math.ceil(value / 10 ** e * 10 ** nd - EPS) / 10 ** nd
        if m >= 10:
            m, e = m / 10, e + 1
        return f"{m:.{nd}f}" == mant and e == int(expo)
    nd = _nd(printed)
    if rule == 'up':
        v = math.ceil(value * 10 ** nd - EPS) / 10 ** nd
    elif rule == 'down':
        v = math.floor(value * 10 ** nd + EPS) / 10 ** nd
    elif rule == 'near':
        v = value
    else:
        raise ValueError(rule)
    return f"{v:.{nd}f}" == printed


def records():
    """(id, key, regex on main.tex + supplement.tex, fn(groups) -> [(printed, value, rule)])."""
    inv = _j('2026-08-28_inverse_moment_falsifier.json')
    dep, summ = inv['deployment'], inv['summary']
    last = dep['sweep'][-1]
    tn = _j('2026-09-27_theory_numerics.json')
    opn = tn['R8_inverse_moment_sensitivity']['open_chain']
    lp = tn['R3_within_sector_lp']['lp_ranges']
    obc = max(c['resid_identity_halfnegT'] / abs(c['m_minus1']) for c in inv['identity_checks'] if not c['pbc'])
    kry = _j('2026-11-08_krylov_detection_order.json')
    taus = [x for r in kry['rows'] for x in r['tau_rel_default_by_k'][1:]]
    pe = _j('2026-08-28_markov_krein_window.json')['orders']['n=4']['peak_edge_window']
    d0 = _j('2026-09-27_delta0_reference_mc.json')
    leak50 = d0['device']['50000']['delta0'] / d0['system']['m0']
    INV = '2026-08-28_inverse_moment_falsifier.json'
    R3 = '2026-09-27_theory_numerics.json:R3_within_sector_lp.lp_ranges'
    P = '(\\d+(?:\\.\\d+)?)'
    return [
        # A5 (main text, Sec. III B and Fig. 3)
        ('A5-text-upto', INV + ':deployment.sweep[-1].{d1,f}',
         r"differs from the reconstruction's by up to \$" + P + r"\\%\$ of \$m_1\$ at \$f=" + P + r"\$",
         lambda g: [(g[0], 100 * last['d1'], 'up'), (g[1], last['f'], 'near')]),
        ('A5-caption-upto', INV + ':deployment.sweep[-1].{d1,d2,dm1}',
         r"move by up to \$" + P + r"\\%\$, \$" + P + r"\\%\$ and \$" + P + r"\\%\$",
         lambda g: [(g[0], 100 * last['d1'], 'up'), (g[1], 100 * last['d2'], 'up'), (g[2], 100 * last['dm1'], 'up')]),
        ('A5-text-compare', INV + ':summary.max_Delta_{minus1,1}_over_m (= sweep[-1].{dm1,d1})',
         r"moves most in relative terms \(up to \$" + P + r"\\%\$ against \$" + P + r"\\%\$ for \$m_1\$\)",
         lambda g: [(g[0], 100 * summ['max_Delta_minus1_over_m'], 'up'), (g[1], 100 * summ['max_Delta_1_over_m'], 'up'),
                    ('summary = sweep[-1]', abs(summ['max_Delta_minus1_over_m'] - last['dm1']) < 1e-12
                     and abs(summ['max_Delta_1_over_m'] - last['d1']) < 1e-12, 'bool')]),
        # A5 (Supplemental Material, Sec. S1 and S8)
        ('A5-obc-identity', INV + ':identity_checks[OBC] + theory_numerics.json:R8.open_chain.resid_poles_vs_fsum_rel',
         r"to a relative \$" + P + r"\\times10\^\{(-\d+)\}\$ at \$L=6\$ and \$8\$ and \$" + P
         + r"\\times10\^\{(-\d+)\}\$ at \$L=12\$",
         lambda g: [((g[0], g[1]), obc, 'sci_up'), ((g[2], g[3]), opn['resid_poles_vs_fsum_rel'], 'sci_up')]),
        ('A5-ring-identity', INV + ':deployment.resid_fsum_vs_poles',
         r"verified at \$L=12\$ to a relative \$" + P + r"\\times10\^\{(-\d+)\}\$",
         lambda g: [((g[0], g[1]), dep['resid_fsum_vs_poles'], 'sci_up')]),
        ('A5-tau-enclosure', '2026-11-08_krylov_detection_order.json:rows[*].tau_rel_default_by_k[1:]',
         r"lies between \$" + P + r"\\%\$ and \$" + P + r"\\%\$ for every \$k=1,\\dots,10\$",
         lambda g: [(g[0], 100 * min(taus), 'down'), (g[1], 100 * max(taus), 'up'),
                    (f'{len(taus)} values', len(taus) == 160, 'bool')]),
        # A1 (Supplemental Material, Sec. S7 and Table SI)
        ('A1-m2-max', R3 + '[0].{Delta_max,max_over_tau}',
         r"the largest shift of \$m_2\$ is \$" + P + r"\$ \(\$" + P + r"\\,\\tau_2\$\)",
         lambda g: [(g[0], lp[0]['Delta_max'], 'near'), (g[1], lp[0]['max_over_tau'], 'near')]),
        ('A1-m3-max', R3 + '[1].{Delta_max,max_over_tau}',
         r"the largest \$m_3\$ shift is \$" + P + r"\$ \(\$" + P + r"\\,\\tau_3\$",
         lambda g: [(g[0], lp[1]['Delta_max'], 'near'), (g[1], lp[1]['max_over_tau'], 'near')]),
        ('A1-tableSI-row1', R3 + '[0].{Delta_min,Delta_max,tau_same_rule,max_over_tau}',
         r"\$m_0,m_1\$ & \$m_2\$ & \$\[(-" + P[1:-1] + r"),\\,"+ P + r"\]\$ & \$" + P + r"\$ & \$" + P + r"\$ &",
         lambda g: [(g[0], lp[0]['Delta_min'], 'down'), (g[1], lp[0]['Delta_max'], 'up'),
                    (g[2], lp[0]['tau_same_rule'], 'near'), (g[3], lp[0]['max_over_tau'], 'near')]),
        ('A1-tableSI-row2', R3 + '[1].{Delta_min,Delta_max,tau_same_rule,max_over_tau}',
         r"\$m_0,\\dots,m_2\$ & \$m_3\$ & \$\[(-" + P[1:-1] + r"),\\," + P + r"\]\$ & \$" + P + r"\$ & \$" + P + r"\$ &",
         lambda g: [(g[0], lp[1]['Delta_min'], 'down'), (g[1], lp[1]['Delta_max'], 'up'),
                    (g[2], lp[1]['tau_same_rule'], 'near'), (g[3], lp[1]['max_over_tau'], 'near')]),
        ('A1-tableSI-row3', R3 + '[2].{Delta_min,Delta_max,tau_same_rule,max_over_tau}',
         r"\$m_0,\\dots,m_3\$ & \$m_4\$ & \$\[(-" + P[1:-1] + r"),\\," + P + r"\]\$ & \$" + P + r"\$ & \$" + P + r"\$ &",
         lambda g: [(g[0], lp[2]['Delta_min'], 'down'), (g[1], lp[2]['Delta_max'], 'up'),
                    (g[2], lp[2]['tau_same_rule'], 'near'), (g[3], max(abs(lp[2]['Delta_min']), lp[2]['Delta_max'])
                                                             / lp[2]['tau_same_rule'], 'near')]),
        # A3 (Supplemental Material, Sec. S6)
        ('A3-mk-edge-window', '2026-08-28_markov_krein_window.json:orders.n=4.peak_edge_window.{interval,width}',
         r"gives \$\\mu\\in\[" + P + r"," + P + r"\]\$, width \$" + P + r"\$",
         lambda g: [(g[0], pe['interval'][0], 'down'), (g[1], pe['interval'][1], 'up'), (g[2], pe['width'], 'near')]),
        # A4 (Supplemental Material, Sec. S10)
        ('A4-leak-50k', '2026-09-27_delta0_reference_mc.json:device.50000.delta0 / system.m0',
         r"carry \$" + P + r"\\%\$ of the weight \(\$1-w=" + P + r"\$\)",
         lambda g: [(g[0], 100 * (1 - leak50), 'near'), (g[1], leak50, 'near')]),
        ('A4-leak-score-50k', '2026-09-27_delta0_reference_mc.json:device.50000.delta0 / system.m0',
         r"\$\\Delta_0/m_0=" + P + r"\$ already flags",
         lambda g: [(g[0], leak50, 'near')]),
        ('A4-error-above', '2026-09-27_delta0_reference_mc.json:device.50000.delta0 / system.m0',
         r"error above \$1-w=" + P + r"\$",
         lambda g: [(g[0], leak50, 'down')]),
    ]


def evaluate(text, recs=None):
    """list of (rid, key, shown, ok); a record whose sentence is absent gives one failing line."""
    out = []
    for rid, key, pat, fn in (recs or records()):
        hits = re.findall(pat, text)
        if not hits:
            out.append((rid, key, 'sentence ABSENT', False))
            continue
        for h in hits:
            g = h if isinstance(h, tuple) else (h,)
            for printed, value, rule in fn(g):
                ok = good(printed, value, rule)
                if rule == 'sci_up':
                    shown = f"{printed[0]}e{printed[1]} vs {value:.4e} ({rule})"
                elif rule == 'bool':
                    shown = f"{printed}: {bool(value)}"
                else:
                    shown = f"{printed} vs {value:.6g} ({rule})"
                out.append((rid, key, shown, ok))
    return out


def failures(text, recs=None):
    return sum(not ok for *_, ok in evaluate(text, recs))


def report(text):
    """print one line per check; return 0 if all pass, 1 otherwise."""
    res = evaluate(text)
    print("\nClosing records (2026-11-14): numbers this text prints in a form the deposit's verify.py does not anchor:")
    for rid, key, shown, ok in res:
        print(f"  [{'PASS' if ok else 'FAIL'}] {rid:<20} {shown}   ({key})")
    nf = sum(not ok for *_, ok in res)
    print(f"  {len(res) - nf} of {len(res)} closing-record values pass")
    print("CLOSING RECORDS: " + ("ALL PASS" if nf == 0 else f"{nf} FAIL"))
    return 1 if nf else 0


if __name__ == '__main__':
    import sys
    PAPER = HERE
    t = ''.join(open(os.path.join(PAPER, f), encoding='utf-8').read() + '\n' for f in ('main.tex', 'supplement.tex'))
    sys.exit(report(t))
