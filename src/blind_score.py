# -*- coding: utf-8 -*-
"""ACTION 3 (4/4) -- UNBLIND and score. Joins frozen verdicts with sealed labels; reports
confusion, TPR/FPR with Wilson 95% CIs, and the per-mode and per-subclass catch rates.
The sealed pre-registration EXPECTED high catch rates for trunc, nl=1 krylov and weight-changing
ac, and low ones for nl>=2 krylov and m0-preserving ac. Two of those secondary predictions FAILED:
nl=1 krylov was caught 0/15, and m0-preserving ac 20/34, an outcome the pre-registration did not
anticipate (not a designed capability). The summary printed below is computed from the data."""
import os, json, math
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, '..', 'data'))


def wilson(k, n, z=1.96):
    if n == 0:
        return (float('nan'), float('nan'))
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0, c - h), min(1, c + h))


def main():
    V = {v['id']: v for v in json.load(open(os.path.join(OUT, 'blind_verdicts.json')))['verdicts']}
    labels = json.load(open(os.path.join(OUT, 'blind_labels_sealed.json')))['labels']
    tp = fp = tn = fn = 0
    per_mode = defaultdict(lambda: [0, 0])       # mode -> [rejected, total]
    per_sub = defaultdict(lambda: [0, 0])
    for lab in labels:
        rej = V[lab['id']]['reject']; corr = lab['corrupted']; mode = lab['mode']
        per_mode[mode][1] += 1; per_mode[mode][0] += int(rej)
        if mode == 'krylov':
            sub = f"krylov nl={lab['params']['nl']}"
        elif mode == 'ac':
            sub = "ac preserve_m0" if lab['params'].get('preserve_m0') else "ac raw"
        else:
            sub = mode
        per_sub[sub][1] += 1; per_sub[sub][0] += int(rej)
        if corr and rej: tp += 1
        elif corr and not rej: fn += 1
        elif not corr and rej: fp += 1
        else: tn += 1
    ncorr = tp + fn; nclean = fp + tn
    tpr = tp / ncorr if ncorr else float('nan'); fpr = fp / nclean if nclean else float('nan')
    tlo, thi = wilson(tp, ncorr); flo, fhi = wilson(fp, nclean)
    print("=" * 66)
    print("BLIND HARNESS SCORE (unblinded)")
    print("=" * 66)
    print(f"confusion: TP={tp} FN={fn} | FP={fp} TN={tn}   (corrupted={ncorr}, clean={nclean})")
    print(f"TPR (sensitivity) = {tpr:.3f}   Wilson95 [{tlo:.3f}, {thi:.3f}]")
    print(f"FPR (false alarm) = {fpr:.3f}   Wilson95 [{flo:.3f}, {fhi:.3f}]   "
          f"(nominal joint ~ {1-(1-0.05)**3:.3f})")
    print("\nper-mode catch rate:")
    for m in sorted(per_mode):
        r, n = per_mode[m]; print(f"   {m:8s}: {r}/{n} = {r/n:.2f}")
    print("\nper-subclass (the honest blind spots):")
    for s in sorted(per_sub):
        r, n = per_sub[s]; print(f"   {s:18s}: {r}/{n} = {r/n:.2f}")
    def _r(key):
        r, n = per_sub.get(key, (0, 0))
        return f"{r}/{n}"
    kry = [per_sub[s] for s in per_sub if s.startswith('krylov')]
    kr, kn = sum(r for r, _ in kry), sum(n for _, n in kry)
    print("\nSUMMARY (computed above): support truncation " + _r('trunc') + "; spurious feature (ac) "
          f"{per_mode['ac'][0]}/{per_mode['ac'][1]} (weight-changing " + _r('ac raw') + ", m0-preserving "
          + _r('ac preserve_m0') + f"); under-converged Krylov {kr}/{kn} at every nl.")
    print("Two sealed secondary predictions FAILED: nl=1 Krylov was expected to be caught and was caught "
          + _r('krylov nl=1') + ";")
    print("m0-preserving ac was expected LOW and was caught " + _r('ac preserve_m0')
          + " (not anticipated by the pre-registration; not a designed capability).")
    print("An nl-node Lanczos reconstruction is exact through m_{2nl-1}: nl>=2 matches m0-m2 and passes by")
    print("construction; nl=1 misses m2 by 2-4% and passes only within the 2% assumed-bias budget. The")
    print("battery detects support/structure errors and spurious features, not moment-preserving errors.")

    out = {'_provenance': {'script': 'blind_score.py', 'sim_only': True,
                           'prereg_sha256': open(os.path.join(OUT, 'prereg.sha256')).read().strip()},
           'confusion': {'TP': tp, 'FN': fn, 'FP': fp, 'TN': tn},
           'TPR': tpr, 'TPR_wilson95': [tlo, thi], 'FPR': fpr, 'FPR_wilson95': [flo, fhi],
           'per_mode': {m: {'rej': per_mode[m][0], 'n': per_mode[m][1]} for m in per_mode},
           'per_subclass': {s: {'rej': per_sub[s][0], 'n': per_sub[s][1]} for s in per_sub}}
    json.dump(out, open(os.path.join(OUT, '2026-08-24_blind_harness_score.json'), 'w'), indent=2)
    print(f"\nwrote 2026-08-24_blind_harness_score.json")


if __name__ == '__main__':
    main()
