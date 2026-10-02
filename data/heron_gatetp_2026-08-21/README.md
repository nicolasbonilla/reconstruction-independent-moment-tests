# `heron_gatetp_2026-08-21/`: raw records of the two-arm ibm_fez run of 2026-08-21

These are the unmodified records of two IBM Quantum SamplerV2 jobs on `ibm_fez`. They were retrieved
read-only from the IBM Quantum Platform REST API (`/jobs/{id}`, `/jobs/{id}/results`, `/jobs/{id}/metrics`)
on 2026-10-02 and copied here byte for byte.

| Job | Arm | U in the R_zz angle | Created (UTC) | Quantum seconds |
|---|---|---|---|---|
| `da44t8u1vhnc73fkie30` | correct | 4.0 | 2026-08-21T13:12:03Z | 69 |
| `da459juaa69c739jcb80` | detuned | 3.5 | 2026-08-21T13:38:23Z | 69 |

Each job has 7 PUBs (K = 7 circuits, k = 0..6), 35 000 shots each, 12 measured bits, and gate and
measurement twirling plus dynamical decoupling switched on. The circuits are those of the
submitting notebook `GATE_TP_on_Heron.ipynb` (cell 4; not deposited): the L = 6 Hubbard ring, an X-gate determinant with 4 up
and 3 down electrons, and then 3 Trotter layers of R_zz(2U dt) and R_xx R_yy hopping, with
dt = k * 0.5 / 3. The job JSON holds the transpiled circuits as they were submitted. When they are
simulated on their 12 active qubits, they reproduce the notebook circuits of their own schedule. This
check is in `src/gate_tp_analysis.py`.

Files per job:
- `<id>_job.json`: the submitted parameters, i.e. the transpiled circuits (QPY, RuntimeEncoder), shots and options.
- `<id>_results.json`: the PrimitiveResult. Each PUB is a BitArray in which each shot is 2 big-endian bytes, and clbit 0 is the least significant bit.
- `<id>_metrics.json`: timestamps and usage.

The files contain no credential: no token, API key, CRN or instance field. The job JSON has a `user_id`
field, which is an IBMid account identifier.

**Bit order.** qiskit writes the counts string with clbit 0 as the rightmost character. The Fock
configuration is therefore `int(bs, 2)`. The notebook read the counts with `int(bs[::-1], 2)`. That
reading maps the (N = 7, S_z = +1/2) sector (4 up, 3 down) onto (N = 7, S_z = -1/2), so every shot it
kept was a device error. (The `Sz` strings in `data/2026-10-02_gatetp_reanalysis.json` give
N_up - N_down, i.e. 2 S_z.)

**Analysis:** `python src/gate_tp_analysis.py` writes `data/2026-10-02_gatetp_reanalysis.json`. No moment
is estimated from these counts: the device supplies only the sampled support.

**Licence:** CC BY 4.0.

SHA-256:
```
a1d8df160b3588d76095f39bd7ea44fae203daee77d4bc65d64be49d39806a73  da44t8u1vhnc73fkie30_job.json
1269b2a574bfbe8b7ac34e17f994f9cfb58192d2e0eb8f1a1dc8a5a22e8474b2  da44t8u1vhnc73fkie30_metrics.json
febf16777b488e5152d089d08837e6d532c046acd131d137fc919b8b2c0f5e44  da44t8u1vhnc73fkie30_results.json
ae148241a88bfead099125f695d0eb7f7363b2d47d7cac539f4817c0164a4950  da459juaa69c739jcb80_job.json
1312e655742f53c63fe4d0e6679bfe0083ba85f5ae5bfd78ec7104465248e0c0  da459juaa69c739jcb80_metrics.json
e91c8cfac11899a13be06446897c322c119d0357a9f3c25f64f052ed60ebcb16  da459juaa69c739jcb80_results.json
```
