# -*- coding: utf-8 -*-
"""Does a manifest-v1-compliant 12-qubit chain exist RIGHT NOW?

Gate v4 denied the amendment, so manifest v1 stands unamended. Its day-of chain rule,
verbatim from manifest_nk_device_v1.json:

    "accept the chain iff min over chain qubits of day-of calibrated T1_q >= 150 us AND
     T2_q >= 100 us AND no chain CZ error > 2x device median; otherwise re-select the
     chain or forfeit the window"

Note what v1 does NOT contain: any readout cap, and any ABSOLUTE CZ cap. Those were
additions proposed by the (denied) amendment, alongside a LOOSENING of T1/T2 to 100/70.
The green-light watcher has been testing the loosened floors this whole time, so it has
never answered the question this script asks.

The canary is already validated at these floors: gate v3R2 returned QUOTABLE-PASS,
10/10, ucb 0.026018 < 0.026667, with the standing caveat that headroom is <= 26%.

Read-only calibration query. ZERO QPU time. Run: python v1_chain_check.py
"""
import sys
import numpy as np

T1_MIN, T2_MIN, NQ = 150e-6, 100e-6, 12
BACKENDS = ('ibm_fez', 'ibm_marrakesh', 'ibm_kingston')


def safe(fn, *a):
    try:
        v = fn(*a)
        return v if v and v > 0 else None
    except Exception:
        return None


def main():
    from qiskit_ibm_runtime import QiskitRuntimeService
    svc = QiskitRuntimeService()

    for name in BACKENDS:
        try:
            be = svc.backend(name)
            props = be.properties()
            cmap = be.coupling_map
        except Exception as e:
            print(f'{name}: no disponible ({type(e).__name__})')
            continue

        n = be.num_qubits
        czs = {}
        for e in cmap.get_edges():
            g = safe(props.gate_error, 'cz', list(e))
            if g is None:
                g = safe(props.gate_error, 'ecr', list(e))
            if g is not None:
                czs[tuple(sorted(e))] = g
        if not czs:
            print(f'{name}: sin datos de puerta de dos qubits')
            continue
        med = float(np.median(list(czs.values())))
        cz_cap = 2.0 * med                      # the ONLY CZ condition in v1

        ok_q = {q for q in range(n)
                if (safe(props.t1, q) or 0) >= T1_MIN and (safe(props.t2, q) or 0) >= T2_MIN}

        adj = {}
        for (a, b), g in czs.items():
            if a in ok_q and b in ok_q and g <= cz_cap:
                adj.setdefault(a, []).append(b)
                adj.setdefault(b, []).append(a)

        # deterministic search for a simple path of NQ qubits, minimising summed CZ error
        best = None

        def walk(path, seen, cost):
            nonlocal best
            if len(path) == NQ:
                if best is None or cost < best[1]:
                    best = (list(path), cost)
                return
            if best is not None and cost >= best[1]:
                return
            for nb in sorted(adj.get(path[-1], []),
                             key=lambda x: czs[tuple(sorted((path[-1], x)))]):
                if nb in seen:
                    continue
                seen.add(nb); path.append(nb)
                walk(path, seen, cost + czs[tuple(sorted((path[-2], path[-1])))])
                path.pop(); seen.discard(nb)

        for s in sorted(ok_q):
            walk([s], {s}, 0.0)
            if best is not None:
                break

        print(f'\n=== {name} ===')
        print(f'  CZ mediana        : {med*1e3:.3f}e-3     tope v1 (2x mediana) = {cz_cap*1e3:.3f}e-3')
        print(f'  qubits con T1>=150us y T2>=100us : {len(ok_q)} de {n}')
        if best is None:
            print(f'  cadena de {NQ} qubits conforme a v1: NO')
        else:
            path, cost = best
            t1s = [safe(props.t1, q) * 1e6 for q in path]
            t2s = [safe(props.t2, q) * 1e6 for q in path]
            ros = [safe(props.readout_error, q) or 0.0 for q in path]
            edges = [czs[tuple(sorted((path[i], path[i+1])))] for i in range(NQ - 1)]
            print(f'  cadena de {NQ} qubits conforme a v1: SI  -> {path}')
            print(f'    T1 min {min(t1s):.0f} us   T2 min {min(t2s):.0f} us')
            print(f'    CZ max en cadena {max(edges)*1e3:.3f}e-3  (tope {cz_cap*1e3:.3f}e-3)')
            print(f'    readout max {max(ros):.4f}   [v1 NO impone tope de readout]')
            print(f'    suma de error CZ {cost*1e3:.3f}e-3')


if __name__ == '__main__':
    main()
