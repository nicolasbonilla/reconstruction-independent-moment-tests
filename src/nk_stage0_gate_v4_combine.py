# -*- coding: utf-8 -*-
"""GATE v4 COMBINER -- applies the sealed amendment rule across the verdict rows.

The rows were executed as separate processes (ROWS=... per process), so the amendment
verdict in nk_stage0_gate_v4.py never fired: it only runs inside the process that owns
row R1 and only when every row is present in that process. This script performs exactly
that cross-row combination, with no new decision rule.

THE SEALED RULE, restated verbatim from nk_stage0_gate_v4.py:
    "Verdict rows (ALL must pass) ... If R1^R2^R3 pass -> manifest v2 reseal;
     if R1 fails solely via readout -> mechanical ladder 0.030->0.025->0.018;
     if R1@0.018 fails or R2 or R3 fail -> amendment DENIED in toto, manifest v1
     stands (gate-negative branch)."
    ONE-AMENDMENT DISCIPLINE: first and LAST floor amendment on this manifest line.
    THR_IV = 0.20*(2/3)/5 = 0.0266667, row passes iff mean + 1.645*SE < THR_IV.

PROVENANCE, stated plainly. Only the final ladder rung has a complete 10-seed row file.
Rungs 0.030 and 0.025 were run in earlier sessions and are recorded here from their run
logs. That does not affect the outcome: the mechanical ladder descends on failure, every
observed replica on both rungs exceeded THR_IV, and the rung that carries the formal
verdict is the last one, which is complete.

HISTORY OF THE SEALED FILE, stated plainly. The first version of this script (2026-09-05) wrote
a `finding` asserting that tightening the readout cap by a factor 1.7 "did not move the canary"
and that the ladder "descended the wrong axis". The ladder in the same file refutes that: the
canary fell 14.2%, against the 18.4% needed to cross the threshold. The same day the `finding`
was rewritten by hand, the caveats on the axis attribution and a `corrections` record (with the
superseded digest) were added, and the file was re-sealed. The committed verdict is therefore the
combiner output plus that hand correction. This version emits the corrected finding, the caveats
and the correction record, with the sealed date, so a recombination in a scratch copy reproduces
the sealed file byte for byte. The corrected finding still names the two-qubit error rate as the
binding constraint; the caveats qualify it (eps is a compound noise knob, not a single physical
parameter).

SIM-ONLY, $0 QPU. Writes data/2026-09-05_nk_v4_verdict.json (+ .sha256); refuses to overwrite
the committed, sealed verdict.
"""
import os, json, glob, hashlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', 'data'))
THR_IV = 0.20 * (2.0 / 3.0) / 5.0                      # 0.0266667, sealed, unchanged
SEALED_DATE = '2026-09-05'                             # date of the sealed verdict (see HISTORY)

# The hand correction of 2026-09-05, recorded in the sealed file (see HISTORY above).
CAVEATS_ON_THE_AXIS_CLAIM = [
    'Only R2 is a clean single-variable comparison against the R1 rung; R3 changes four '
    'parameters at once and cannot carry the attribution on its own.',
    "eps is a compound knob in this harness: theta_zz = sqrt(eps) moves the depolarizing and the "
    "coherent-ZZ components together, so 'the two-qubit axis' is not a single physical parameter.",
    'the readout axis is applied through an oracle-exact inversion, so the two axes are not on '
    'equal footing as noise models.',
]
CORRECTIONS = [{
    'date': '2026-09-05',
    'what': 'the `finding` field asserted that tightening readout did not move the canary; the '
            'ladder in this same file refutes it (-14.19% achieved vs -18.41% needed). Rewritten, '
            'and the caveats on the axis attribution added.',
    'and': 'the original .sha256 was computed over an LF stream while the file was written CRLF, '
           'so it failed verification. Re-sealed over the bytes on disk.',
    'superseded_digest': 'f1b7aa486fa50f03b3b29a9a6d3ed609355347813ec3c34873f4e0396838be3d',
}]

# ladder rungs 1 and 2, recovered from their run logs (see PROVENANCE above)
LADDER_PARTIAL = {
    'R1(ro=0.030)': [0.0355, 0.0301, 0.0337, 0.0331, 0.0340, 0.0297],
    'R1ro0.025':    [0.0322, 0.0302, 0.0319, 0.0298],
}


def load_rows():
    rows = {}
    for f in sorted(glob.glob(os.path.join(RES, '*nk_v4_row_*.json'))):
        tag = os.path.basename(f).split('nk_v4_row_')[1].rsplit('.json', 1)[0]
        rows[tag] = json.load(open(f))
    return rows


def main():
    rows = load_rows()
    print('=' * 74)
    print('GATE v4 -- CROSS-ROW COMBINATION (sealed rule, no new criteria)')
    print('=' * 74)
    print(f'sealed iv-D threshold THR_IV = {THR_IV:.7f}\n')

    print('--- mechanical readout ladder on R1 (eps = 3.5e-3) ---')
    ladder = []
    for tag, vals in LADDER_PARTIAL.items():
        a = np.array(vals)
        ucb = float(a.mean() + 1.645 * a.std(ddof=1) / np.sqrt(len(a)))
        dec = ucb < THR_IV
        ladder.append({'rung': tag, 'n_replicas': len(vals), 'complete': False,
                       'ivD_mean': float(a.mean()), 'ivD_ucb': ucb, 'row_pass': bool(dec),
                       'provenance': 'run log, earlier session (partial)'})
        print(f'  {tag:<16} n={len(vals):>2} (partial)  mean={a.mean():.5f}  '
              f'ucb={ucb:.5f}  -> {"PASS" if dec else "FAIL"}   '
              f'[min replica {a.min():.4f}, all above THR: {bool((a > THR_IV).all())}]')

    r1f = rows.get('R1ro0.018')
    assert r1f is not None, 'final ladder rung R1ro0.018 missing'
    ladder.append({'rung': 'R1ro0.018', 'n_replicas': len(r1f['replicas']), 'complete': True,
                   'ivD_mean': r1f['ivD_mean'], 'ivD_ucb': r1f['ivD_ucb'],
                   'row_pass': r1f['row_pass'], 'provenance': 'complete row file'})
    print(f'  {"R1ro0.018":<16} n={len(r1f["replicas"]):>2} (COMPLETE) mean={r1f["ivD_mean"]:.5f}  '
          f'ucb={r1f["ivD_ucb"]:.5f}  -> {"PASS" if r1f["row_pass"] else "FAIL"}   '
          f'fails={r1f["fail_channels"]}')
    print('  => the sealed ladder 0.030 -> 0.025 -> 0.018 is EXHAUSTED.\n')

    print('--- the other two verdict rows ---')
    for tag in ('R2', 'R3'):
        r = rows[tag]
        p = r['params']
        print(f'  {tag:<16} eps={p["eps"]:.5g} ro={p["ro"]:.3g}  '
              f'ucb={r["ivD_ucb"]:.5f} -> {"PASS" if r["row_pass"] else "FAIL"}')

    amendment_pass = bool(r1f['row_pass'] and rows['R2']['row_pass'] and rows['R3']['row_pass'])
    verdict = ('AMENDMENT-VALIDATED' if amendment_pass
               else 'AMENDMENT-DENIED (manifest v1 stands; gate-negative)')

    # what the ladder establishes, computed from it (the corrected finding of 2026-09-05; see
    # HISTORY in the module docstring): readout moved the canary, but not far enough
    top, bot = ladder[0]['ivD_mean'], r1f['ivD_mean']
    achieved = (top - bot) / top                       # 0.1419
    needed = (top - THR_IV) / top                      # 0.1841
    finding = (
        'The readout ladder is exhausted without passing, but not because readout is inert. '
        f'Tightening the per-qubit readout cap from 0.030 to 0.018 moved the iv-D canary from '
        f'{top:.6f} to {bot:.6f}, a reduction of {100*achieved:.1f}% against the {100*needed:.1f}% '
        f'needed to cross the sealed threshold {THR_IV:.7f}: readout delivered '
        f'{100*achieved/needed:.0f}% of the required movement and then reached its own floor. '
        'The residual is available on the two-qubit axis, and both rows that pass sit there '
        '(R2 at eps = 2.85e-3, R3 at eps = 2.0e-3) rather than at any readout value. So the '
        'ladder was not wasted; it was simply one axis short, and the binding constraint for a '
        'future window is the two-qubit error rate.'
    )

    frac = [r['ivD_mean'] / THR_IV for r in (r1f, rows['R2'], rows['R3'])]
    headroom = (f'Across the completed verdict rows the iv-D canary runs at '
                f'{min(frac)*100:.0f}-{max(frac)*100:.0f}% of its sealed threshold.')

    print('\n' + '=' * 74)
    print(f'GATE v4 VERDICT: {verdict}')
    print('=' * 74)
    print(headroom)
    print('\n' + finding)
    print('\nCaveats on the axis attribution:')
    for c in CAVEATS_ON_THE_AXIS_CLAIM:
        print('  - ' + c)

    out = {
        '_provenance': {
            'script': 'nk_stage0_gate_v4_combine.py',
            'date': SEALED_DATE,
            'caveats_on_the_axis_claim': CAVEATS_ON_THE_AXIS_CLAIM,
            'corrections': CORRECTIONS,
            'sim_only': True, 'qpu_spent': 0,
            'sealed_rule': 'R1^R2^R3 must pass; ladder 0.030->0.025->0.018 on readout-only '
                           'failure; exhausted ladder => amendment denied in toto',
            'one_amendment_discipline': 'first and LAST floor amendment on this manifest '
                                        'line; a second is forbidden threshold-shopping',
            'partial_rungs_note': 'rungs 0.030 and 0.025 are recorded from run logs of '
                                  'earlier sessions and are not complete 10-seed rows; every '
                                  'observed replica on both exceeded THR_IV, and the rung '
                                  'carrying the formal verdict (0.018) is complete',
        },
        'THR_IV': THR_IV,
        'ladder': ladder,
        'rows': {t: {'params': rows[t]['params'], 'ivD_mean': rows[t]['ivD_mean'],
                     'ivD_ucb': rows[t]['ivD_ucb'], 'row_pass': rows[t]['row_pass'],
                     'fail_channels': rows[t]['fail_channels']} for t in sorted(rows)},
        'amendment_pass': amendment_pass,
        'verdict': verdict,
        'headroom_line': headroom,
        'finding': finding,
    }
    path = os.path.join(RES, '2026-09-05_nk_v4_verdict.json')
    if os.path.exists(path) or os.path.exists(path + '.sha256'):   # sealed, read-only
        raise SystemExit('REFUSING to overwrite the sealed verdict ' + os.path.basename(path)
                         + '; run in a scratch copy of the repository to recombine.')
    blob = json.dumps(out, indent=1, sort_keys=True)
    open(path, 'w', newline='\n').write(blob)   # LF on every OS, so the file hashes to its seal
    digest = hashlib.sha256(blob.encode()).hexdigest()
    open(path + '.sha256', 'w', newline='\n').write(digest + '  ' + os.path.basename(path) + '\n')
    print(f'\nwrote {os.path.basename(path)}')
    print(f'sha256 {digest}')


if __name__ == '__main__':
    main()
