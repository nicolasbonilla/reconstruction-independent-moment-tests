# -*- coding: utf-8 -*-
"""How far is the hardware from a manifest-v1-compliant 12-qubit chain?

v1_chain_check.py answers yes/no. That is not actionable on its own: 'no' could mean the
longest qualifying chain is 11 (one qubit away, and a matter of waiting for a
recalibration) or 4 (structurally out of reach). This reports the LONGEST qualifying
chain on each backend, and attributes the binding constraint to T1/T2 or to the CZ
homogeneity condition, so the answer is a number rather than a verdict.

v1 rule: T1 >= 150 us, T2 >= 100 us per chain qubit; no chain CZ > 2x device median.
No readout cap and no absolute CZ cap exist in v1.

Read-only calibration query. ZERO QPU. Run: python v1_chain_diagnose.py
"""
import numpy as np

T1_MIN, T2_MIN, NQ = 150e-6, 100e-6, 12
BACKENDS = ('ibm_kingston', 'ibm_fez', 'ibm_marrakesh')


def safe(fn, *a):
    try:
        v = fn(*a)
        return v if v and v > 0 else None
    except Exception:
        return None


def longest_path(adj, cap=NQ):
    """Longest simple path, capped at `cap`; DFS from every node with memo-free pruning."""
    best = [0, None]

    def walk(u, seen, path):
        if len(path) > best[0]:
            best[0], best[1] = len(path), list(path)
        if len(path) >= cap:
            return
        for v in adj.get(u, []):
            if v not in seen:
                seen.add(v); path.append(v)
                walk(v, seen, path)
                path.pop(); seen.discard(v)

    for s in sorted(adj):
        walk(s, {s}, [s])
        if best[0] >= cap:
            break
    return best


def main():
    from qiskit_ibm_runtime import QiskitRuntimeService
    svc = QiskitRuntimeService()

    for name in BACKENDS:
        be = svc.backend(name)
        props, cmap = be.properties(), be.coupling_map
        n = be.num_qubits

        czs = {}
        for e in cmap.get_edges():
            g = safe(props.gate_error, 'cz', list(e)) or safe(props.gate_error, 'ecr', list(e))
            if g is not None:
                czs[tuple(sorted(e))] = g
        med = float(np.median(list(czs.values())))
        cap = 2.0 * med

        okq = {q for q in range(n)
               if (safe(props.t1, q) or 0) >= T1_MIN and (safe(props.t2, q) or 0) >= T2_MIN}

        # three graphs: T1/T2 only, CZ only, and both (the v1 rule)
        def build(qfilter, efilter):
            a = {}
            for (x, y), g in czs.items():
                if qfilter(x) and qfilter(y) and efilter(g):
                    a.setdefault(x, []).append(y)
                    a.setdefault(y, []).append(x)
            return a

        g_t = build(lambda q: q in okq, lambda g: True)          # T1/T2 only
        g_c = build(lambda q: True, lambda g: g <= cap)          # CZ homogeneity only
        g_b = build(lambda q: q in okq, lambda g: g <= cap)      # the v1 rule

        Lt = longest_path(g_t)[0]
        Lc = longest_path(g_c)[0]
        Lb, pb = longest_path(g_b)

        print(f'\n=== {name} ===')
        print(f'  CZ mediana {med*1e3:.3f}e-3   tope v1 = 2x mediana = {cap*1e3:.3f}e-3')
        print(f'  qubits que pasan T1>=150us y T2>=100us : {len(okq)}/{n}')
        print(f'  aristas que pasan el tope de CZ        : '
              f'{sum(1 for g in czs.values() if g <= cap)}/{len(czs)}')
        print(f'  cadena mas larga -- solo T1/T2         : {Lt}')
        print(f'  cadena mas larga -- solo homogeneidad CZ: {Lc}')
        print(f'  cadena mas larga -- REGLA v1 COMPLETA   : {Lb}   '
              f'{"(SUFICIENTE)" if Lb >= NQ else f"(faltan {NQ-Lb})"}')
        if Lb >= NQ:
            print(f'    cadena: {pb[:NQ]}')
        else:
            binding = 'T1/T2' if Lt < Lc else ('homogeneidad CZ' if Lc < Lt else 'ambas')
            print(f'    restriccion vinculante: {binding}')


if __name__ == '__main__':
    main()
