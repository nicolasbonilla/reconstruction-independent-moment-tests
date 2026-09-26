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

SIM-ONLY, $0 QPU. Writes data/2026-09-05_nk_v4_verdict.json (+ .sha256); refuses to overwrite
the committed, sealed verdict.
"""
import os, json, glob, hashlib, datetime
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', 'data'))
THR_IV = 0.20 * (2.0 / 3.0) / 5.0                      # 0.0266667, sealed, unchanged

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

    # what the ladder actually establishes: the binding axis is the two-qubit error, not readout
    finding = (
        'The readout ladder is exhausted without passing: at eps = 3.5e-3 the iv-D canary '
        'exceeds its sealed threshold at ro = 0.030, 0.025 and 0.018 alike, while both rows '
        'that pass require a LOWER two-qubit error (R2 at eps = 2.85e-3, R3 at eps = 2.0e-3). '
        'The binding constraint is therefore the two-qubit error rate, not the readout cap: '
        'tightening readout by a factor 1.7 did not move the canary, and lowering eps by a '
        'factor 1.2 crossed it. The pre-registered ladder descended the wrong axis.'
    )

    frac = [r['ivD_mean'] / THR_IV for r in (r1f, rows['R2'], rows['R3'])]
    headroom = (f'Across the completed verdict rows the iv-D canary runs at '
                f'{min(frac)*100:.0f}-{max(frac)*100:.0f}% of its sealed threshold.')

    print('\n' + '=' * 74)
    print(f'GATE v4 VERDICT: {verdict}')
    print('=' * 74)
    print(headroom)
    print('\n' + finding)

    out = {
        '_provenance': {
            'script': 'nk_stage0_gate_v4_combine.py',
            'date': datetime.date.today().isoformat(),
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
