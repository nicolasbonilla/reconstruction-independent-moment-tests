# -*- coding: utf-8 -*-
"""Green-light watcher (0 QPU): every CHECK_S seconds test (a) usage_remaining >= 120s and
(b) the sealed day-of chain rule finds a 12-qubit chain on the CURRENT fez calibration.
Exits GREEN the moment both hold (the sealed job is then fired by the operator);
exits AMBER after MAX_CHECKS without a green light."""
import sys, time, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
CHECK_S, MAX_CHECKS = 5400, 16
from qiskit_ibm_runtime import QiskitRuntimeService
from nk_stage2_device_job import pick_chain
for i in range(MAX_CHECKS):
    try:
        svc = QiskitRuntimeService()
        u = svc.usage()
        rem = u['usage_remaining_seconds']
        backend = svc.backend('ibm_fez')
        cal = backend.properties().last_update_date
        best, med = pick_chain(backend)
        chain_ok = best is not None
        print(f"[check {i+1}/{MAX_CHECKS}] {time.strftime('%H:%M')} remaining={rem}s cal={cal} chain={'YES' if chain_ok else 'no'}", flush=True)
        if chain_ok and rem >= 120:
            print(f"GREEN: chain={best[0]} mean_cz={best[1]:.4f} remaining={rem}s", flush=True)
            sys.exit(0)
    except Exception as e:
        print(f"[check {i+1}] error: {str(e)[:150]}", flush=True)
    time.sleep(CHECK_S)
print("AMBER: no green light within the watch window", flush=True)
