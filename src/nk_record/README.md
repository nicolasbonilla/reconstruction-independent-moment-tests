# n_k device attempt record (2026-09-03 to 2026-09-05)

The scripts of the sealed on-device n_k falsifier attempt, as they were run, plus the generator of
`data/2026-09-05_nk_device_attempt_record.json`. No device job was ever submitted under manifest v1
(0 jobs, 0 shots); the record states why, with the logged outputs and a re-query of the calibrations.

| file | what it is |
|---|---|
| `nk_stage2_device_job_launch1.py` | the job as launched first on 2026-09-03 (byte-identical to commit 9bdb1d5) |
| `nk_stage2_device_job_launch2.py` | the job as launched second (launch 1 plus two edits, rebuilt from the session log) |
| `nk_stage2_device_job_20260903_1505Z.py` | the module the watchers imported (exhaustive search, same v1 rule) |
| `green_light_watch.py`, `_v2.py` | calibration watchers for the v1 rule on ibm_fez (observe only) |
| `green_light_watch_v3.py` | watcher for the proposed amended rule (never sealed; observe only) |
| `v1_chain_check.py`, `v1_chain_diagnose.py` | the 2026-09-05 checks of the sealed v1 rule on three backends |
| `nk_logged_outputs.py` | the logged outputs, transcribed verbatim with their UTC timestamps |
| `nk_attempt_record.py` | builds the record; `--offline` recomputes it from the stored calibration tables |

These historical scripts keep their original paths (`../06_results`) and imports; they are a record, not
part of the pipeline. Run the generator with `python src/nk_record/nk_attempt_record.py --offline`.
The data record is released under CC BY 4.0.
