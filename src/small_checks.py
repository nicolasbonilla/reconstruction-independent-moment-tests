# -*- coding: utf-8 -*-
r"""small_checks.py -- small committed theory checks for the P3 revision, and the single writer of
data/2026-09-27_theory_numerics.json (one key per plan item; each key carries its own provenance).

Keys written by THIS script (the others are written by the scripts listed in KEY_OWNERS):
  R11_trace_over_dim   Tr(O)/dim for O = c_{0,up} (H-E0) c_{0,up}^dag on the L=8, U/t=4 periodic chain of the
                       device job (Hamiltonian rebuilt exactly as in hardware_matched_job_L8.py), computed as a
                       numerical trace and compared with the closed form U(1/4+(L-1)/8) - E0/2; m1 = <0|O|0>;
                       trace-to-value ratio; first-order global-depolarizing bias p[Tr(O)/dim - m1].
  m4_psd_and_symmetry  (i) is O positive semidefinite for the bare-H convention? (lowest eigenvalue of O);
                       (ii) critic item 12: m0 = 1-<n_{0,up}> = 1/2 because <n_{j,up}> = N_up/L by translation
                       invariance (checked on the half-filled ground state and, at another filling, on the
                       ground-manifold average), not specifically by particle-hole symmetry;
                       (iii) main.tex:453 '<rho_pi> vanishes by half-filling symmetry' tested on the ground state
                       that bond_moment_estimator.py actually uses (open chain, global Fock-space ground state).
  M5_christoffel_claims, M6_inverse_moment_claim   derived summaries of keys R4/R9 and R8 (no new physics).

Run the whole A3 set (writes every key; 14-15 min on the authors' laptop in the 2026-09-27 runs, of which
christoffel_tolerance_lp.py, export_christoffel_dat.py, necessary_sufficient_composition.py,
run_inverse_moment_falsifier.py --r8 and this script's own R11/m4 keys take 2-3 min each):
                                                                 cd src && python small_checks.py --all
Run only this script's keys:                                     cd src && python small_checks.py
No random numbers are drawn here (seed = none); eigensolvers are dense or started from a fixed vector.
"""
import os, sys, json, time, datetime, platform, subprocess

import numpy as np
import scipy
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = os.path.normpath(os.path.join(HERE, '..', 'data'))
THEORY_JSON = os.path.join(DATA, '2026-09-27_theory_numerics.json')

KEY_OWNERS = {
    'R3_within_sector_lp': 'within_sector_lp.py',
    'M1_shifted_power_bound': 'within_sector_lp.py',
    'R4_tolerance_inflated_atom': 'christoffel_tolerance_lp.py',
    'R5_hausdorff_weighted_support': 'necessary_sufficient_composition.py',
    'R6_estimator_forms_off_eigenstate': 'estimator_form_offeigenstate.py',
    'R8_inverse_moment_sensitivity': 'run_inverse_moment_falsifier.py --r8',
    'R9_christoffel_max_Wn': 'export_christoffel_dat.py',
    'm1_radau_outside_hull': 'export_christoffel_dat.py',
    'm2_rescaled_frame': 'export_christoffel_dat.py',
    'R11_trace_over_dim': 'small_checks.py',
    'm4_psd_and_symmetry': 'small_checks.py',
    'M5_christoffel_claims': 'small_checks.py (derived from R4, R9)',
    'M6_inverse_moment_claim': 'small_checks.py (derived from R8)',
}
RUN_ORDER = ['within_sector_lp.py', 'christoffel_tolerance_lp.py', 'necessary_sufficient_composition.py',
             'estimator_form_offeigenstate.py', 'run_inverse_moment_falsifier.py --r8',
             'export_christoffel_dat.py']


# ------------------------------------------------------------------------------------------------
# shared writer
# ------------------------------------------------------------------------------------------------
def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return [_jsonable(v) for v in o.tolist()]
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, complex):
        return {'re': o.real, 'im': o.imag}
    return o


def _env():
    import mpmath
    return {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
            'mpmath': mpmath.__version__, 'platform': platform.platform()}


def record_theory_numerics(key, payload, script, runtime_s, seed=None, command=None):
    """Insert/replace one key of data/2026-09-27_theory_numerics.json (read-modify-atomic-write).
    The file name is fixed (it is the dated record of this revision step); the run date is stored inside."""
    doc = {}
    if os.path.exists(THEORY_JSON):
        with open(THEORY_JSON, 'r', encoding='utf-8') as fh:
            doc = json.load(fh)
    payload = dict(payload)
    payload['provenance'] = {'script': 'src/' + script, 'command': command or f'cd src && python {script}',
                             'seed': seed, 'runtime_s': round(float(runtime_s), 2),
                             'run_date': datetime.date.today().isoformat(), 'environment': _env()}
    doc[key] = payload
    doc['provenance'] = {
        'file': 'data/2026-09-27_theory_numerics.json',
        'what': 'P3 revision (SciPost Physics Core), area A3 theory numerics: plan items R3 R4 R5 R6 R8 R9 R11 '
                'M1 M5 M6 m1 m2 m4 and critic items 5, 12 (PLAN_P3_2026-09-26.md). One key per item; each key '
                'has its own provenance (script, seed, runtime). Simulation/exact numerics only; no device data.',
        'regenerate_all': 'cd src && python small_checks.py --all',
        'key_owners': KEY_OWNERS,
        'keys_present': sorted(k for k in doc if k != 'provenance'),
        'last_write': datetime.datetime.now().isoformat(timespec='seconds'),
    }
    tmp = THEORY_JSON + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(_jsonable(doc), fh, indent=1, ensure_ascii=False)
    os.replace(tmp, THEORY_JSON)
    print(f"[theory_numerics] wrote key '{key}' -> data/2026-09-27_theory_numerics.json")


def load_theory_numerics():
    with open(THEORY_JSON, 'r', encoding='utf-8') as fh:
        return json.load(fh)


# ------------------------------------------------------------------------------------------------
# R11 + m4(i),(ii): the L=8, U/t=4 device system (Hamiltonian as in hardware_matched_job_L8.py)
# ------------------------------------------------------------------------------------------------
def _popcount(x, nbits):
    c = np.zeros_like(x)
    for b in range(nbits):
        c += (x >> b) & 1
    return c


def fock_L8(L=8, U=4.0, thop=1.0):
    """Fock-space c-operators and H exactly as hardware_matched_job_L8.py:61-83 (interleaved p=2*site+spin,
    Jordan-Wigner sign over all lower modes, periodic chain, bare U n_up n_dn). Vectorized, same matrices."""
    M = 2 * L; dim = 1 << M
    s = np.arange(dim, dtype=np.int64)
    C = []
    for p in range(M):
        occ = ((s >> p) & 1).astype(bool)
        sign = (-1.0) ** _popcount(s[occ] & ((1 << p) - 1), M)
        C.append(sp.csr_matrix((sign, (s[occ] & ~(1 << p), s[occ])), shape=(dim, dim)))
    Cd = [c.T.conj().tocsr() for c in C]
    Nocc = _popcount(s, M)
    Szocc = np.zeros(dim, dtype=np.int64)
    for i in range(L):
        Szocc += ((s >> (2 * i)) & 1) - ((s >> (2 * i + 1)) & 1)
    H = sp.csr_matrix((dim, dim))
    for i in range(L):
        j = (i + 1) % L
        for spin in (0, 1):
            a = 2 * i + spin; b = 2 * j + spin
            H = H - thop * (Cd[a] @ C[b] + Cd[b] @ C[a])
    for i in range(L):
        H = H + U * (Cd[2 * i] @ C[2 * i]) @ (Cd[2 * i + 1] @ C[2 * i + 1])
    return C, Cd, H.tocsr(), Nocc, Szocc, dim


def check_R11_m4_device_system():
    t0 = time.time()
    L, U = 8, 4.0
    C, Cd, H, Nocc, Szocc, dim = fock_L8(L, U)
    # half-filled ground state (Nocc==L, Szocc==0) -- as the device job
    gi = np.where((Nocc == L) & (Szocc == 0))[0]
    eg, vg = np.linalg.eigh(H[gi][:, gi].toarray())
    E0 = float(eg[0]); gap_sector = float(eg[1] - eg[0])
    psi0 = np.zeros(dim); psi0[gi] = vg[:, 0]
    phi = Cd[0] @ psi0
    m0 = float(np.vdot(phi, phi).real)
    m1 = float(np.vdot(phi, H @ phi).real - E0 * m0)
    # --- R11: Tr(O)/dim as a numerical trace ---
    O = (C[0] @ (H - E0 * sp.identity(dim, format='csr')) @ Cd[0]).tocsr()
    tr_over_dim = float(O.diagonal().sum().real) / dim
    closed = U * (0.25 + (L - 1) / 8.0) - E0 / 2.0
    # committed device record for m1 (read-only cross-check)
    dev = json.load(open(os.path.join(DATA, 'heron_counts_matched_L8_da983qse74ec73ajfr2g.json')))
    m1_committed = float(dev['verification']['m1_op'])
    # --- m4(ii): translation invariance of <n_{j,sigma}> (critic 12) ---
    n_modes = np.array([float(np.vdot(psi0, (Cd[p] @ C[p]) @ psi0).real) for p in range(2 * L)])
    n_up = n_modes[0::2]; n_dn = n_modes[1::2]
    # another filling: (N_up,N_dn)=(3,3) -> ground-manifold-averaged density
    other = {}
    for (nu, nd) in [(3, 3), (5, 3)]:
        Nn, Sz = nu + nd, nu - nd
        ii = np.where((Nocc == Nn) & (Szocc == Sz))[0]
        e, v = np.linalg.eigh(H[ii][:, ii].toarray())
        deg = int(np.sum(e < e[0] + 1e-8))
        dens_single = []; dens_avg = np.zeros(2 * L)
        for g in range(deg):
            vec = np.zeros(dim); vec[ii] = v[:, g]
            dn = np.array([float(np.vdot(vec, (Cd[p] @ C[p]) @ vec).real) for p in range(2 * L)])
            dens_avg += dn / deg
            if g == 0:
                dens_single = dn
        other[f'Nup{nu}_Ndn{nd}'] = {
            'ground_degeneracy': deg, 'E_ground': float(e[0]),
            'manifold_avg_n_up_per_site': dens_avg[0::2], 'manifold_avg_n_dn_per_site': dens_avg[1::2],
            'expected_N_up_over_L': nu / L, 'expected_N_dn_over_L': nd / L,
            'max_dev_manifold_avg_from_Nsigma_over_L': float(max(np.max(np.abs(dens_avg[0::2] - nu / L)),
                                                                 np.max(np.abs(dens_avg[1::2] - nd / L)))),
            'single_eigh_vector_n_up_per_site': dens_single[0::2],
            'single_vector_translation_invariant': bool(np.max(np.abs(dens_single[0::2] - nu / L)) < 1e-8),
            'addition_m0_at_site0_up_manifold_avg': float(1.0 - dens_avg[0]),
        }
    # --- m4(i): is O PSD? lowest eigenvalue of the compression of H onto 'mode 0 occupied', minus E0 ---
    occ0 = ((np.arange(dim) >> 0) & 1).astype(bool)
    worst = (np.inf, None); per_sector = []
    for Nn in range(1, 2 * L + 1):
        for Sz in range(-L, L + 1):
            ii = np.where(occ0 & (Nocc == Nn) & (Szocc == Sz))[0]
            if len(ii) == 0:
                continue
            Hb = H[ii][:, ii]
            if len(ii) <= 400:
                e = float(np.linalg.eigvalsh(Hb.toarray())[0])
            else:
                v0 = np.random.default_rng(0).standard_normal(len(ii))
                e = float(eigsh(Hb, k=1, which='SA', v0=v0, tol=1e-10)[0][0])
            per_sector.append((Nn, Sz, len(ii), e))
            if e < worst[0]:
                worst = (e, (Nn, Sz, len(ii)))
    lam_min_O = min(0.0, worst[0] - E0)
    runtime = time.time() - t0
    R11 = {
        'system': 'L=8, U/t=4 periodic chain, bare U n_up n_dn, interleaved JW order (hardware_matched_job_L8.py); '
                  'O = c_{0,up}(H-E0)c_{0,up}^dag, E0 = half-filled (N_up,N_dn)=(4,4) ground energy; dim = 2^16',
        'E0': E0, 'Tr_O_over_dim_numerical': tr_over_dim, 'Tr_O_over_dim_closed_form': closed,
        'closed_form': 'U*(1/4 + (L-1)/8) - E0/2  (Tr(c H c^dag)=Tr(H n_0up): hopping traceless; '
                       'U-term 1/4 on site 0, 1/8 on each other site; Tr(c c^dag)/dim = 1/2)',
        'abs_diff_numerical_vs_closed': abs(tr_over_dim - closed),
        'm0_addition': m0, 'm1_addition': m1, 'm1_committed_device_record': m1_committed,
        'abs_diff_m1_vs_committed': abs(m1 - m1_committed),
        'trace_to_value_ratio': tr_over_dim / m1,
        'global_depolarizing_first_order': {
            'model': '<O>_p = (1-p) m1 + p Tr(O)/dim  =>  bias = p [Tr(O)/dim - m1]',
            'bias_per_unit_p_absolute': tr_over_dim - m1,
            'bias_per_unit_p_relative_to_m1': (tr_over_dim - m1) / m1,
            'sign': 'upward (Tr(O)/dim > m1); this does NOT require O to be PSD (see m4_psd_and_symmetry)'},
        'plan_expectation': 'Tr(O)/dim = 6.80 with E0 = -4.6035; main.tex:189 quotes ~6.8, m1 ~2.1, ratio ~3.2',
    }
    m4 = {
        'i_psd_of_O': {
            'question': "main.tex:189 'O = c(H-E0)c^dag is positive semidefinite' for the bare-H convention",
            'lowest_eig_of_H_compressed_to_mode0_occupied': worst[0], 'attained_in_sector_N_Sz_dim': worst[1],
            'E0_half_filled': E0, 'lambda_min_O': lam_min_O,
            'O_is_PSD': bool(lam_min_O >= -1e-9),
            'n_blocks_scanned': len(per_sector),
            'note': 'lambda_min(O) = min(0, min_{N,Sz} lambda_min(P0 H P0) - E0), P0 = projector on mode (0,up) '
                    'occupied. With bare U n_up n_dn and no chemical potential the half-filled E0 is not the global '
                    'Fock minimum, so O has negative eigenvalues; the upward depolarizing bias follows from '
                    'Tr(O)/dim > m1 alone (R11).'},
        'ii_m0_half_reason': {
            'critic_item': 12,
            'n_up_per_site_halffilled_GS': n_up, 'n_dn_per_site_halffilled_GS': n_dn,
            'N_up_over_L': 4 / L,
            'max_dev_from_N_up_over_L': float(max(np.max(np.abs(n_up - 0.5)), np.max(np.abs(n_dn - 0.5)))),
            'ground_state_sector_gap': gap_sector,
            'm0 = 1 - <n_0up>': m0,
            'other_fillings': other,
            'conclusion': 'm0 = 1 - <n_{0,up}> = 1 - N_up/L: the reason is translation invariance of the ground '
                          '(manifold) density, valid at any filling (e.g. 5/8 at (3,3)); particle-hole symmetry is '
                          'not needed. At half filling N_up/L = 1/2.'},
    }
    return R11, m4, runtime


# ------------------------------------------------------------------------------------------------
# m4(iii): <rho_pi> on the estimator's actual ground state (bond_moment_estimator.py)
# ------------------------------------------------------------------------------------------------
def check_m4_rho_shift():
    t0 = time.time()
    os.environ.setdefault('BME_L', '6')
    out = {}
    from scipy.sparse.linalg import eigsh as _eigsh
    import importlib
    for Lx in (4, 6, 8):
        os.environ['BME_L'] = str(Lx)
        import bond_moment_estimator as bme
        bme = importlib.reload(bme)                         # L is read from BME_L at import
        Hm = bme.build_hubbard().to_matrix(sparse=True).tocsr()
        rho = bme.build_rho().to_matrix(sparse=True).tocsr()
        NQ = 2 * Lx; dim = 1 << NQ
        s = np.arange(dim, dtype=np.int64)
        nq = [((s >> q) & 1).astype(float) for q in range(NQ)]      # qubit q occupied (blocked order)
        Ndiag = np.sum(nq, axis=0)
        v0 = np.random.default_rng(0).standard_normal(dim)
        ev, V = _eigsh(Hm, k=6, which='SA', v0=v0, tol=1e-12)
        o = np.argsort(ev); ev, V = ev[o], V[:, o]
        deg = int(np.sum(ev < ev[0] + 1e-8))
        # the state build_operators() uses: eigsh k=1 (take our lowest vector; for a degenerate ground
        # manifold we report every vector and the manifold average)
        rows = []
        for g in range(deg):
            psi = V[:, g]; p = np.abs(psi) ** 2
            nsite = np.array([float(np.sum(p * (nq[i] + nq[Lx + i]))) for i in range(Lx)])
            rows.append({'N': float(np.sum(p * Ndiag)), 'rho_pi': float(np.vdot(psi, rho @ psi).real),
                         'n_site': nsite, 'reflection_max_dev': float(np.max(np.abs(nsite - nsite[::-1])))})
        # half-filled sector ground state, for contrast
        hf = np.where(np.abs(Ndiag - Lx) < 1e-9)[0]
        Hhf = Hm[hf][:, hf]
        eh, vh = _eigsh(Hhf, k=1, which='SA', v0=np.random.default_rng(1).standard_normal(len(hf)), tol=1e-12)
        psih = np.zeros(dim, dtype=complex); psih[hf] = vh[:, 0]
        out[f'L{Lx}'] = {
            'global_ground_energy': float(ev[0]), 'global_ground_degeneracy': deg,
            'global_ground_vectors': rows,
            'global_ground_is_half_filled': bool(all(abs(r['N'] - Lx) < 1e-6 for r in rows)),
            'half_filled_sector_E0': float(eh[0]),
            'half_filled_sector_rho_pi': float(np.vdot(psih, rho @ psih).real),
        }
    os.environ['BME_L'] = '6'
    runtime = time.time() - t0
    return {
        'question': "main.tex:453 '<rho_q> vanishes by half-filling symmetry' (also estimator_scaling.py docstring)",
        'estimator_setup': 'bond_moment_estimator.py: OBC chain, blocked JW order, U n_up n_dn written in Paulis, '
                           'q = pi; build_operators() takes the GLOBAL Fock-space ground state (eigsh, no chemical '
                           'potential, no particle-number constraint)',
        'by_L': out,
        'mechanism': 'for even L the bond-centred reflection j -> L-1-j maps cos(pi j) -> -cos(pi j); a '
                     'reflection-symmetric density therefore gives <rho_pi> = 0 at ANY filling',
    }, runtime


# ------------------------------------------------------------------------------------------------
# derived summaries
# ------------------------------------------------------------------------------------------------
def summaries():
    doc = load_theory_numerics()
    out = {}
    if 'R4_tolerance_inflated_atom' in doc and 'R9_christoffel_max_Wn' in doc:
        R4 = doc['R4_tolerance_inflated_atom']; R9 = doc['R9_christoffel_max_Wn']
        out['M5_christoffel_claims'] = {
            'derived_from': ['R4_tolerance_inflated_atom', 'R9_christoffel_max_Wn', 'm2_rescaled_frame'],
            'deployed_order': 'n=1 (m0,m1,m2)',
            'max_t_W1_exact': R9['exact_max']['1']['max_W'],
            'max_t_W1_equals_m0_at_mean': R9['W1_at_mean_equals_m0'],
            'max_t_W4_exact': R9['exact_max']['4']['max_W'],
            'max_t_Wn_exact': {n: R9['exact_max'][n]['max_W'] for n in R9['exact_max']},
            'n_distinct_atoms_with_weight': R9['n_weighted_atoms'],
            'W_n_at_dominant_atom': R9['W_at_dominant_atom'],
            'dominant_residue': R9['dominant_residue'],
            'W_n_at_dominant_atom_equals_residue_for_n_ge_Nminus1': bool(
                abs(R9['W_at_dominant_atom'][str(R9['n_weighted_atoms'] - 1)] - R9['dominant_residue']) < 1e-9),
            'U4_rows_t_12_15_20_inflated_atom': {r['label']: r['inflated_atom_closed_R']
                                                 for r in R4['U4_within_sector']['rows']
                                                 if r['label'] in ('t=12', 't=15', 't=20')},
            'tolerance_inflated_over_W1_U4': {r['label']: r['ratio_inflated_over_W1']
                                              for r in R4['U4_within_sector']['rows']},
            'tolerance_inflated_over_W1_U8': {r['label']: r['ratio_inflated_over_W1']
                                              for r in R4['U8_christoffel_measure']['rows_rule_tau']},
            'claim_20_30_moments': 'unsupported: the measure has ' + str(R9['n_weighted_atoms']) +
                                   ' distinct weighted atoms; W_n at the dominant atom saturates at its residue for '
                                   'n >= N-1 and never goes below it, so no number of moments gives a ~0.1 bracket '
                                   'at the peak',
        }
    if 'R8_inverse_moment_sensitivity' in doc:
        R8 = doc['R8_inverse_moment_sensitivity']
        ring = R8['ring']; obc = R8['open_chain']
        out['M6_inverse_moment_claim'] = {
            'derived_from': ['R8_inverse_moment_sensitivity'],
            'relative_shift_at_f0.5_ring': {'m_-1': ring['rel_shift_at_f'][-1]['dm1'],
                                            'm1': ring['rel_shift_at_f'][-1]['d1'],
                                            'm2': ring['rel_shift_at_f'][-1]['d2']},
            'ratio_m-1_over_m1_relative_shift_ring': ring['slope_ratio_minus1_over_1'],
            'ratio_m-1_over_m2_relative_shift_ring': ring['slope_ratio_minus1_over_2'],
            'ring_cancellation_D_over_half_negT': ring['cancellation_D_over_half_negT'],
            'f_star_ring_bias_only': ring['f_star_bias_only'],
            'f_star_ring_full_rule': ring['f_star_full_rule'],
            'f_star_open_chain_bias_only': obc['f_star_bias_only'],
            'f_star_open_chain_full_rule': obc['f_star_full_rule'],
            'relative_slope_open_chain': obc['slope_relative_per_unit_f'],
        }
        verdict = {}
        for tag, rec in (('ring', ring), ('open_chain', obc)):
            rows = [(rec['f_star_full_rule'], 2.5)] + [(s['f_star_full_rule'], s['om_spur'])
                                                       for s in rec.get('om_spur_scan', [])]
            wins = [osp for fs, osp in rows if fs['m_-1'] is not None and fs['m1'] is not None and fs['m2'] is not None
                    and fs['m_-1'] < min(fs['m1'], fs['m2'])]
            verdict[tag] = {'om_spur_where_m_-1_fires_first_full_rule': sorted(set(wins)),
                            'om_spur_scanned': sorted(set(osp for _, osp in rows))}
            if 'om_spur_scan' in rec:
                out['M6_inverse_moment_claim'][f'om_spur_scan_{tag}_f_star_full_rule'] = {
                    str(s['om_spur']): s['f_star_full_rule'] for s in rec['om_spur_scan']}
        out['M6_inverse_moment_claim']['noise_normalized_verdict'] = verdict
    return out


# ------------------------------------------------------------------------------------------------
def main():
    if '--all' in sys.argv:
        for cmd in RUN_ORDER:
            print(f"\n===== python {cmd} =====", flush=True)
            r = subprocess.run([sys.executable] + cmd.split(), cwd=HERE)
            if r.returncode != 0:
                raise SystemExit(f"{cmd} failed with code {r.returncode}")
    print("===== small_checks: R11 + m4 (L=8 device system) =====", flush=True)
    R11, m4, rt = check_R11_m4_device_system()
    print(f"Tr(O)/dim numerical = {R11['Tr_O_over_dim_numerical']:.6f}  closed form = "
          f"{R11['Tr_O_over_dim_closed_form']:.6f}  E0 = {R11['E0']:.6f}  m1 = {R11['m1_addition']:.6f}  "
          f"ratio = {R11['trace_to_value_ratio']:.3f}")
    print(f"lambda_min(O) = {m4['i_psd_of_O']['lambda_min_O']:.6f}  (PSD: {m4['i_psd_of_O']['O_is_PSD']}); "
          f"<n_j,up> max dev from N_up/L = {m4['ii_m0_half_reason']['max_dev_from_N_up_over_L']:.2e}")
    record_theory_numerics('R11_trace_over_dim', R11, 'small_checks.py', rt, seed=None)
    rho, rt2 = check_m4_rho_shift()
    for k, v in rho['by_L'].items():
        print(f"{k}: global GS deg={v['global_ground_degeneracy']} N={[round(r['N'],6) for r in v['global_ground_vectors']]} "
              f"rho_pi={[round(r['rho_pi'],12) for r in v['global_ground_vectors']]}  half-filled sector rho_pi="
              f"{v['half_filled_sector_rho_pi']:.3e}")
    m4['iii_rho_pi_shift'] = rho
    record_theory_numerics('m4_psd_and_symmetry', m4, 'small_checks.py', rt + rt2, seed=0,
                           command='cd src && python small_checks.py  (eigsh start vectors: default_rng(0/1))')
    t0 = time.time()
    for k, v in summaries().items():
        record_theory_numerics(k, v, 'small_checks.py', time.time() - t0)


if __name__ == '__main__':
    main()
