# -*- coding: utf-8 -*-
"""Green-light watcher v3 — implements day_of_rule_v2 (amended, wf_c07c723b-b96) as OBSERVER.
Backends B={fez, marrakesh, kingston}; qualification CZmed<=3.5e-3; chain rule: 12-path with
T1>=100us, T2>=70us, RO<=RO_CAP per qubit, edge CZ<=min(3.5e-3, 2*med). Polls calibration
timestamps every 15 min; full check on any NEW calibration. Exits GREEN when a chain exists
AND budget >=130s. Reports both RO caps 0.030 and 0.025 (ladder-pending; firing uses the
SEALED value only, after manifest v2). 0 QPU."""
import sys, time, os, collections
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from qiskit_ibm_runtime import QiskitRuntimeService
B = ('ibm_fez', 'ibm_marrakesh', 'ibm_kingston')
def safe(fn,*a):
    try: return fn(*a)
    except Exception: return None
def chain(bk, props, ro_cap):
    czd = {}
    for (x,y) in bk.coupling_map:
        e = safe(props.gate_error,'cz',[x,y])
        if e is not None: czd[(x,y)] = e
    med = float(np.median(list(czd.values())))
    if med > 3.5e-3: return None, med
    cap = min(3.5e-3, 2*med)
    ok = {q for q in range(bk.num_qubits)
          if (safe(props.t1,q) or 0)>=100e-6 and (safe(props.t2,q) or 0)>=70e-6
          and (safe(props.readout_error,q) or 1)<=ro_cap}
    adj = collections.defaultdict(set)
    for (x,y),e in czd.items():
        if x in ok and y in ok and e<=cap: adj[x].add(y); adj[y].add(x)
    def dfs(n,vis,path):
        if len(path)>=12: return list(path)
        for m in sorted(adj[n]):
            if m not in vis:
                vis.add(m); path.append(m)
                r = dfs(m,vis,path)
                if r: return r
                vis.remove(m); path.pop()
        return None
    for s in sorted(ok):
        r = dfs(s,{s},[s])
        if r: return r, med
    return None, med
svc = QiskitRuntimeService()
last = {}
for i in range(144):
    try:
        found_any = False
        for bn in B:
            bk = svc.backend(bn)
            props = bk.properties(refresh=True)
            cal = str(props.last_update_date)
            sealed_now = os.path.exists(os.path.join(HERE, '..', '06_results',
                                                     'manifest_nk_device_v2.json'))
            if cal == last.get(bn) and not sealed_now: continue
            last[bn] = cal
            c30, med = chain(bk, props, 0.030)
            c25, _   = chain(bk, props, 0.025)
            rem = svc.usage()['usage_remaining_seconds']
            print(f"[{time.strftime('%d %H:%M')}] {bn} NEW cal {cal[-8:]} CZmed={med*1e3:.2f} "
                  f"chain@030={'YES' if c30 else 'no'} chain@025={'YES' if c25 else 'no'} rem={rem}s", flush=True)
            if c30 and rem >= 130:
                sealed = os.path.exists(os.path.join(HERE, '..', '06_results',
                                                     'manifest_nk_device_v2.json'))
                if sealed:
                    print(f"GREEN(ro030): {bn} chain={c30} rem={rem}s", flush=True)
                    found_any = True
                else:
                    print(f"  [green window OPEN on {bn} but manifest v2 not sealed -> "
                          f"cannot fire; chain={c30}]", flush=True)
        if found_any: sys.exit(0)
        if i % 8 == 0: print(f"[{time.strftime('%d %H:%M')}] pulse {i+1}/144 (no new/green)", flush=True)
    except Exception as e:
        print(f"[pulse {i+1}] err {str(e)[:120]}", flush=True)
        try: svc = QiskitRuntimeService()
        except Exception: pass
    time.sleep(900)
print("AMBER: 36h without green", flush=True)
