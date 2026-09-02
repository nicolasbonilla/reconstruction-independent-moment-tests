# -*- coding: utf-8 -*-
"""make_paper_figs.py — all 6 paper figures in the ART-FIRST publication style (paper_style.py).
Reuses the verified result JSONs in 06_results; recomputes only what is cheap. Beautiful, colorful, cohesive."""
import os, sys, json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.normpath(os.path.join(HERE, '..', '..', '00_shared', 'lib'))
sys.path.insert(0, LIB); sys.path.insert(0, HERE)
import paper_style as PS; PS.apply()
RES = os.path.normpath(os.path.join(HERE, '..', '06_results'))
PAP = os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs')); os.makedirs(PAP, exist_ok=True)
FIGSRC = os.path.normpath('C:/Users/Nicolas/Downloads/Proyecto_SQD_ML/02_Paper_Amortizacion/release/paper/figs')


def load(name):
    return json.load(open(os.path.join(RES, name)))


def save(fig, stem):
    fig.savefig(os.path.join(PAP, stem + '.pdf')); fig.savefig(os.path.join(PAP, stem + '.png'), dpi=200)
    plt.close(fig); print('saved', stem)


# ---------- FIG 1 : hero schematic ----------
def fig_hero():
    fig, ax = plt.subplots(figsize=(7.0, 3.95)); ax.set_xlim(0, 10); ax.set_ylim(-0.7, 6); ax.axis('off')
    def box(x, y, w, h, fc, ec, txt, fs=11, tc=PS.INK):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.06,rounding_size=0.12',
                                    fc=fc, ec=ec, lw=1.4, mutation_aspect=1))
        ax.text(x + w / 2, y + h / 2, txt, ha='center', va='center', fontsize=fs, color=tc)
    def arrow(x1, y1, x2, y2, col):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=15,
                                     lw=1.8, color=col, shrinkA=2, shrinkB=2))
    box(0.3, 2.4, 2.2, 1.2, '#EAF2F8', PS.SLATE, 'computational-\nbasis samples\n$\\mathcal{S}$', 10.5)
    # top branch: reconstruction
    box(3.6, 4.1, 2.9, 1.3, '#FDEDE7', PS.WARMDEEP, 'reconstruction\n$A(\\omega)$', 11)
    # bottom branch: independent moments
    box(3.6, 0.6, 2.9, 1.3, '#E8F5F1', PS.TEAL,
        'independent moments\n$\\widehat m_k=\\langle O_k\\rangle$', 10.5)
    box(7.6, 2.4, 2.1, 1.2, '#FBEFC7', PS.WARM, 'compare\n$\\Delta_k=|\\widehat m_k-\\bar m_k|$', 9.8)
    arrow(2.5, 3.2, 3.6, 4.6, PS.WARMDEEP)
    arrow(2.5, 2.8, 3.6, 1.2, PS.TEAL)
    arrow(6.5, 4.5, 7.7, 3.3, PS.WARMDEEP)
    arrow(6.5, 1.2, 7.7, 2.7, PS.TEAL)
    ax.text(8.65, 1.75, 'falsify  $\\Delta_k\\!>\\!\\epsilon$', ha='center', fontsize=9.5, color=PS.CRIMSON, style='italic')
    ax.text(8.65, 4.35, 'corroborate\n$\\Delta_k\\!\\approx\\!0$', ha='center', fontsize=9.5, color=PS.TEAL, style='italic')
    arrow(8.6, 3.55, 8.6, 4.0, PS.TEAL); arrow(8.6, 2.4, 8.6, 2.05, PS.CRIMSON)
    ax.text(5.0, 5.75, 'the reconstruction-independent moment loop', ha='center', fontsize=13, color=PS.INK)
    ax.text(5.0, -0.45, r'with $\ m_0=\langle J^2\rangle,\ \ m_1=\frac{1}{2}\langle[J,[H,J]]\rangle,\ \ O_k=k$-fold commutator', ha='center',
            fontsize=9, color=PS.MUTE, style='italic')
    save(fig, 'fig_hero')


# ---------- FIG 3 : teeth ----------
def fig_teeth():
    d = load('2026-08-18_falsifier_teeth.json')['sweep']
    cov = [100 * r['coverage'] for r in d]
    ri = [100 * r['residual_independent'] for r in d]
    rc = [100 * r['residual_circular'] for r in d]
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.fill_between(cov, 0, ri, color=PS.WARM, alpha=0.14)
    ax.plot(cov, ri, 'o-', color=PS.WARMDEEP, mfc=PS.WARM, mec=PS.WARMDEEP, mew=1.1,
            label='independent estimator $\\Delta_k$ — fires')
    ax.plot(cov, rc, 's--', color=PS.TEAL, mfc='white', mec=PS.TEAL, mew=1.3,
            label='back-computed control — blind')
    ax.axhline(5, color=PS.CRIMSON, ls=':', lw=1.3)
    ax.text(96, 8, 'falsification threshold', fontsize=8.5, color=PS.CRIMSON, style='italic')
    ax.invert_xaxis()
    ax.set_xlabel('spectral coverage kept  (\\%)'.replace('\\%', '%'))
    ax.set_ylabel('first-moment discrepancy  (\\%)'.replace('\\%', '%'))
    ax.set_title('Independence gives the screen teeth', fontsize=12.5)
    ax.legend(loc='upper left'); PS.finish(ax); save(fig, 'fig_teeth')


# ---------- FIG 4 : gauss law ----------
def fig_gausslaw():
    d = load('2026-08-18_gausslaw_falsifier.json')['error_sweep']
    p = [100 * r['p'] for r in d]
    g = [r['gauss_violation'] for r in d]
    q = [r['charge_drift'] for r in d]
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.fill_between(p, 0, g, color=PS.TEAL, alpha=0.13)
    ax.plot(p, g, 'o-', color=PS.TEAL, mfc='#8FD3C4', mec=PS.TEAL, mew=1.1,
            label=r'Gauss-law falsifier $\langle\sum_x G_x^2\rangle$ — fires')
    ax.plot(p, q, 's--', color=PS.SLATE, mfc='white', mec=PS.SLATE, mew=1.3,
            label='global charge check — blind')
    ax.set_xlabel('per-qubit gauge-breaking error rate  (\\%)'.replace('\\%', '%'))
    ax.set_ylabel('violation signal')
    ax.set_title('Cross-domain: an exact Gauss-law falsifier', fontsize=12.5)
    ax.legend(loc='upper left'); PS.finish(ax); save(fig, 'fig_gausslaw')


# ---------- FIG 5 : Heron ----------
def fig_heron():
    d = load('2026-08-18_heron_screen.json')
    sweep = d['teeth_on_real_data_distortion_sweep']
    ex = np.loadtxt(os.path.join(FIGSRC, 'heron_hot.dat'), skiprows=1)
    hw = np.loadtxt(os.path.join(FIGSRC, 'heron_hw.dat'), skiprows=1)
    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.4)); fig.subplots_adjust(wspace=0.28)
    ax[0].fill_between(ex[:, 0], 0, ex[:, 1], color=PS.SLATE, alpha=0.16)
    ax[0].plot(ex[:, 0], ex[:, 1], '-', color=PS.SLATE, lw=1.9, label='exact reconstruction (target)')
    ax[0].plot(hw[:, 0], hw[:, 1], 'o-', color=PS.WARMDEEP, mfc=PS.WARM, mec=PS.WARMDEEP, mew=0.8, ms=4.5, lw=1.2,
               label='real IBM Heron $A(\\omega)$')
    ax[0].set_xlim(2, 12); ax[0].set_xlabel(r'frequency  $\omega-E_0$'); ax[0].set_ylabel(r'$A(\omega)$')
    ax[0].set_title('The screen on real Heron data', fontsize=12); ax[0].legend(loc='upper right'); PS.finish(ax[0])
    f = [100 * s['f'] for s in sweep]
    ax[1].fill_between(f, 0, [100 * s['r1_firstmoment'] for s in sweep], color=PS.WARM, alpha=0.14)
    ax[1].plot(f, [100 * s['r1_firstmoment'] for s in sweep], 'o-', color=PS.WARMDEEP, mfc=PS.WARM, mec=PS.WARMDEEP, mew=1.1,
               label='first moment — fires')
    ax[1].plot(f, [100 * s['r0_weight'] for s in sweep], 's--', color=PS.TEAL, mfc='white', mec=PS.TEAL, mew=1.3,
               label='$f$-sum — blind')
    ax[1].axhline(5, color=PS.CRIMSON, ls=':', lw=1.3)
    ax[1].set_xlabel('induced weight misplacement  (\\%)'.replace('\\%', '%'))
    ax[1].set_ylabel('sum-rule residual  (\\%)'.replace('\\%', '%'))
    ax[1].set_title('Teeth on real hardware data', fontsize=12); ax[1].legend(loc='upper left'); PS.finish(ax[1])
    save(fig, 'fig_heron')


# ---------- FIG 6 : christoffel ----------
def fig_christoffel():
    import hubbard_ed as H, scipy.sparse as sp
    from scipy.sparse.linalg import eigsh
    from run_sumrule_falsifier import current_operator
    from run_moment_bound_theorem import christoffel_widths
    L, U = 6, 8.0
    c = H.build_operators(L); cd = [ci.getH() for ci in c]
    Ham0 = H.hubbard_hamiltonian(L, t=1.0, U=U, pbc=True, c=c)
    Ntot = sum((cd[m] @ c[m] for m in range(2 * L)), sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
    J = current_operator(L, c)
    for mu in [0.0, -1.0, -2.0]:
        cand = eigsh(Ham0 - mu * Ntot, k=1, which='SA')[1][:, 0]
        if float(np.real(np.vdot(cand, Ntot @ cand))) <= L - 1.5:
            psi0 = cand; break
    E0 = float(np.real(np.vdot(psi0, Ham0 @ psi0)))
    Jpsi = J @ psi0
    Ev, Vv = np.linalg.eigh(Ham0.toarray()); om = Ev - E0; w = np.abs(Vv.conj().T @ Jpsi) ** 2
    keep = om > 1e-6; om, w = om[keep], w[keep]; wn = w / w.sum()
    a, b = om.min(), om.max(); xs = (2 * om - (a + b)) / (b - a)
    xq = np.linspace(-0.999, 0.999, 60)
    Wc = christoffel_widths(xs, wn, [1, 2, 3, 4], xq)
    Wc = {n: np.clip(W, 0, 1) for n, W in Wc.items() if np.all(np.isfinite(W)) and float(np.nanmax(W)) > 1e-6}
    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.4)); fig.subplots_adjust(wspace=0.30)
    cols = [PS.SLATE, PS.TEAL, PS.WARM, PS.WARMDEEP]
    for col, (n, W) in zip(cols, sorted(Wc.items())):
        ax[0].fill_between(xq, 0, W, color=col, alpha=0.10)
        ax[0].plot(xq, W, lw=2.0, color=col, label=f'$n={n}$  (moments $\\leq{2*n}$)')
    ax[0].set_xlabel(r'rescaled frequency  $t\in[-1,1]$')
    ax[0].set_ylabel(r'certified miss-distance  $W_n(t)$')
    ax[0].set_title(r'Christoffel band  $W_n(t)=1/\sum_{j\leq n}p_j(t)^2$', fontsize=11.5)
    ax[0].legend(loc='upper right', fontsize=8.5); PS.finish(ax[0])
    ns = sorted(Wc); wm = [float(np.nanmax(Wc[n])) for n in ns]
    ax[1].plot(ns, wm, 'o-', color=PS.WARMDEEP, mfc=PS.WARM, mec=PS.WARMDEEP, mew=1.2, ms=9)
    ax[1].set_xlabel('polynomial order  $n$'); ax[1].set_ylabel(r'max band width  $\max_t W_n(t)$')
    ax[1].set_title('More moments $\\Rightarrow$ tighter certificate', fontsize=11.5)
    ax[1].set_xticks(ns); PS.finish(ax[1])
    save(fig, 'fig_christoffel')


# ---------- FIG 2 : falsifier (recompute) ----------
def fig_falsifier():
    import hubbard_ed as H, scipy.sparse as sp
    from scipy.sparse.linalg import eigsh
    from run_sumrule_falsifier import current_operator
    L, U = 6, 8.0
    c = H.build_operators(L); cd = [ci.getH() for ci in c]
    Ham0 = H.hubbard_hamiltonian(L, t=1.0, U=U, pbc=True, c=c)
    Ntot = sum((cd[m] @ c[m] for m in range(2 * L)), sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
    J = current_operator(L, c)
    for mu in [0.0, -1.0, -2.0]:
        cand = eigsh(Ham0 - mu * Ntot, k=1, which='SA')[1][:, 0]
        if float(np.real(np.vdot(cand, Ntot @ cand))) <= L - 1.5:
            psi0 = cand; break
    E0 = float(np.real(np.vdot(psi0, Ham0 @ psi0)))
    Jpsi = J @ psi0
    m1_op = float(np.real(np.vdot(Jpsi, Ham0 @ Jpsi)) - E0 * np.real(np.vdot(Jpsi, Jpsi)))
    Ev, Vv = np.linalg.eigh(Ham0.toarray()); om = Ev - E0; w = np.abs(Vv.conj().T @ Jpsi) ** 2
    keep = om > 1e-6; om, w = om[keep], w[keep]
    om_split = 4.0; lo, hi = om < om_split, om >= om_split; W_lo, W_hi = w[lo].sum(), w[hi].sum()
    fs = np.linspace(0.0, 0.5, 11); r0, r1 = [], []
    for f in fs:
        ww = w.copy(); ww[hi] *= (1 - f); ww[lo] *= (1 + f * W_hi / max(W_lo, 1e-12))
        r0.append(abs(ww.sum() - w.sum()) / w.sum()); r1.append(abs((om * ww).sum() - m1_op) / m1_op)
    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.5)); fig.subplots_adjust(wspace=0.28)
    bins = np.linspace(0, om.max() * 1.02, 34); fw = 0.5
    ww = w.copy(); ww[hi] *= (1 - fw); ww[lo] *= (1 + fw * W_hi / max(W_lo, 1e-12))
    ax[0].hist(om, bins=bins, weights=w, color=PS.SLATE, alpha=0.85, edgecolor='white', linewidth=0.4, label=r'true $A_J(\omega)$')
    ax[0].hist(om, bins=bins, weights=ww, color=PS.WARMDEEP, alpha=0.60, edgecolor='white', linewidth=0.4, label='wrong split (same $m_0$)')
    ax[0].axvline(om_split, color=PS.INK, ls=(0, (2, 2)), lw=1.0)
    ax[0].set_title(r'A wrong split that \emph{passes} the $f$-sum'.replace('\\emph{', '').replace('}', ''), fontsize=12)
    ax[0].set_xlabel(r'$\omega$  (energy above ground state)'); ax[0].set_ylabel(r'spectral weight $A_J(\omega)$')
    ax[0].legend(loc='upper right'); PS.finish(ax[0])
    ax[1].fill_between(fs, 0, [100 * x for x in r1], color=PS.WARM, alpha=0.16)
    ax[1].plot(fs, [100 * x for x in r1], 'o-', color=PS.WARMDEEP, mfc=PS.WARM, mec=PS.WARMDEEP, mew=1.1, label=r'first moment $m_1$ — fires')
    ax[1].plot(fs, [100 * x for x in r0], 's--', color=PS.TEAL, mfc='white', mec=PS.TEAL, mew=1.3, label=r'$f$-sum $m_0$ — blind')
    ax[1].axhline(5, color=PS.CRIMSON, ls=':', lw=1.3)
    ax[1].set_title(r'The first moment \emph{falsifies} what the $f$-sum misses'.replace('\\emph{', '').replace('}', ''), fontsize=11.5)
    ax[1].set_xlabel(r'weight fraction moved  mid-IR $\rightarrow$ Drude'); ax[1].set_ylabel('sum-rule residual  (\\%)'.replace('\\%', '%'))
    ax[1].legend(loc='upper left'); PS.finish(ax[1])
    save(fig, 'fig_falsifier')


if __name__ == '__main__':
    fig_hero(); fig_falsifier(); fig_teeth(); fig_gausslaw(); fig_heron(); fig_christoffel()
    print('ALL FIGURES DONE ->', PAP)
