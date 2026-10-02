# -*- coding: utf-8 -*-
"""Green-light watcher v2 (0 QPU): every 15 min poll ONLY the fez calibration timestamp
(one cheap call). When a NEW calibration lands, run the full sealed chain-rule DFS +
usage check. Exit GREEN the moment both pass; exit AMBER after ~36 h without green.
One full check also runs at startup (in case conditions changed since the last look)."""
import sys, time, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
PULSE_S, MAX_PULSES = 900, 144            # 15 min x 144 = 36 h
from qiskit_ibm_runtime import QiskitRuntimeService
from nk_stage2_device_job import pick_chain

def full_check(svc):
    u = svc.usage(); rem = u['usage_remaining_seconds']
    backend = svc.backend('ibm_fez')
    best, med = pick_chain(backend)
    return rem, best

svc = QiskitRuntimeService()
last_cal = None
for i in range(MAX_PULSES):
    try:
        backend = svc.backend('ibm_fez')
        cal = str(backend.properties(refresh=True).last_update_date)
        fresh = (cal != last_cal)
        if fresh:
            rem, best = full_check(svc)
            print(f"[{time.strftime('%d %H:%M')}] NEW calibration {cal} -> chain={'YES' if best else 'no'} remaining={rem}s", flush=True)
            if best is not None and rem >= 120:
                print(f"GREEN: chain={best[0]} mean_cz={best[1]:.4f} remaining={rem}s cal={cal}", flush=True)
                sys.exit(0)
            last_cal = cal
        elif i % 8 == 0:
            print(f"[{time.strftime('%d %H:%M')}] pulse {i+1}/{MAX_PULSES}: same calibration ({cal})", flush=True)
    except Exception as e:
        print(f"[pulse {i+1}] error: {str(e)[:150]}", flush=True)
        try: svc = QiskitRuntimeService()
        except Exception: pass
    time.sleep(PULSE_S)
print("AMBER: no green light within the 36 h watch window", flush=True)
