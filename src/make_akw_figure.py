# -*- coding: utf-8 -*-
"""make_akw_figure.py — the RICH showpiece figure: the single-particle spectral function A(k,omega) of the
Hubbard ring (Mott gap, upper/lower Hubbard bands, spinon-holon dispersion), and the moment screen
operating on it. Beautiful warm heatmaps in the paper style, echoing the companion paper's A(k,omega)."""
import os, sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.normpath(os.path.join(HERE, '..', '..', '00_shared', 'lib'))
sys.path.insert(0, LIB); sys.path.insert(0, HERE)
import hubbard_ed as H
import paper_style as PS; PS.apply()
PAP = os.path.normpath(os.path.join(HERE, '..', 'paper', 'figs')); os.makedirs(PAP, exist_ok=True)

L, U = 6, 8.0
ETA = 0.28                       # Lorentzian broadening
NW = 400


def ck(kj, spin, c, L):
    """momentum-space annihilation c_{k,spin} = (1/sqrt L) sum_j e^{-i k j} c_{j,spin}."""
    k = 2 * np.pi * kj / L
    op = None
    for j in range(L):
        term = np.exp(-1j * k * j) * c[2 * j + spin]
        op = term if op is None else op + term
    return op / np.sqrt(L)


def akw():
    c = H.build_operators(L); cd = [ci.getH() for ci in c]
    Ham = H.hubbard_hamiltonian(L, t=1.0, U=U, pbc=True, c=c)
    N = sum((cd[m] @ c[m] for m in range(2 * L)), sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
    Sz = sum((0.5 * (cd[2 * s] @ c[2 * s] - cd[2 * s + 1] @ c[2 * s + 1]) for s in range(L)),
             sp.csr_matrix((2 ** (2 * L),) * 2, dtype=complex))
    I = sp.identity(2 ** (2 * L), dtype=complex, format='csr')
    # half-filling N=L, Sz=0 ground state
    Hpen = Ham + 30 * ((N - L * I) @ (N - L * I)) + 30 * (Sz @ Sz)
    E0, V0 = eigsh(Hpen, k=1, which='SA'); psi0 = V0[:, 0]
    E0 = float(np.real(np.vdot(psi0, Ham @ psi0)))
    Ev, Vv = np.linalg.eigh(Ham.toarray())
    ks = list(range(L + 1))                      # k index 0..L (fold to [0,2] pi via extension)
    wgrid = np.linspace(-9, 9, NW)
    A = np.zeros((len(ks), NW))
    for ik, kj in enumerate(ks):
        ckop = ck(kj % L, 0, c, L)               # spin up
        rem = Vv.conj().T @ (ckop @ psi0)        # <n| c_k |0>  (removal, omega<0)
        add = Vv.conj().T @ (ckop.getH() @ psi0) # <n| c_k^dag |0>  (addition, omega>0)
        wr = np.abs(rem) ** 2; wa = np.abs(add) ** 2
        omr = -(Ev - E0); oma = (Ev - E0)
        for w0, wt in ((omr, wr), (oma, wa)):
            for e, weight in zip(w0, wt):
                if weight > 1e-9:
                    A[ik] += weight * (ETA / np.pi) / ((wgrid - e) ** 2 + ETA ** 2)
    return np.array(ks) / L * 2, wgrid, A     # k in units of pi (0..2)


def main():
    kk, ww, A = akw()
    # a "wrong reconstruction": blur + shift weight across the Mott gap (an SQD-truncation-like distortion)
    from scipy.ndimage import gaussian_filter1d
    Awrong = A.copy()
    Awrong = gaussian_filter1d(Awrong, 6, axis=1)
    # move a slice of upper-band weight down into the gap (wrong)
    mid = NW // 2
    Awrong[:, mid - 30:mid] += 0.35 * Awrong[:, mid + 10:mid + 40]
    Awrong[:, mid + 10:mid + 40] *= 0.65

    # per-k first-moment residual of true vs wrong (the screen)
    def m1(M):
        w = M / (M.sum(axis=1, keepdims=True) + 1e-12)
        return (w * ww[None, :]).sum(axis=1)
    r = np.abs(m1(Awrong) - m1(A))

    fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.6),
                           gridspec_kw={'width_ratios': [1, 1, 0.58], 'wspace': 0.34})
    vmax = np.percentile(A, 99.7)
    im = None
    for col, (Ad, ttl) in enumerate([(A, r'true $A(k,\omega)$'), (Awrong, 'wrong reconstruction')]):
        im = ax[col].imshow(Ad.T, origin='lower', aspect='auto', cmap=PS.WARMMAP, vmin=0, vmax=vmax,
                            extent=[kk.min(), kk.max(), ww.min(), ww.max()])
        ax[col].axhline(0, color='white', lw=1.0, ls=(0, (4, 3)), alpha=0.85)
        ax[col].set_title(ttl, fontsize=12, color=PS.INK)
        ax[col].set_xlabel(r'momentum  $k/\pi$')
        if col == 0:
            ax[col].set_ylabel(r'frequency  $\omega/t$')
            ax[col].text(0.06, 6.7, 'upper Hubbard band', color='white', fontsize=8.5, style='italic')
            ax[col].text(0.06, -7.7, 'lower Hubbard band', color='white', fontsize=8.5, style='italic')
            ax[col].annotate('Mott gap', xy=(1.0, 0), xytext=(1.12, 2.7), color='white', fontsize=9,
                             arrowprops=dict(arrowstyle='<->', color='white', lw=1.1))
        else:
            ax[col].set_yticklabels([])
    cb = fig.colorbar(im, ax=ax[1], fraction=0.05, pad=0.03)
    cb.set_label(r'$A(k,\omega)$', fontsize=10); cb.ax.tick_params(labelsize=8)
    ax[2].fill_betweenx(kk, 0, 100 * r, color=PS.WARM, alpha=0.16)
    ax[2].plot(100 * r, kk, 'o-', color=PS.WARMDEEP, mfc=PS.WARM, mec=PS.WARMDEEP, mew=1.0, ms=5)
    ax[2].axvline(5, color=PS.CRIMSON, ls=':', lw=1.2)
    ax[2].set_title('the screen\nfires per $k$', fontsize=10.5)
    ax[2].set_xlabel(r'first-moment residual (\%)'.replace('\\%', '%'))
    ax[2].yaxis.set_label_position('right'); ax[2].yaxis.tick_right(); ax[2].set_ylabel(r'$k/\pi$')
    ax[2].grid(alpha=0.2); PS.finish(ax[2])
    ax[2].spines['right'].set_visible(True); ax[2].spines['left'].set_visible(False)
    fig.suptitle(r'The moment screen on the single-particle spectral function  $A(k,\omega)$  of the Hubbard model',
                 y=1.00, fontsize=11.5)
    fig.tight_layout()
    fig.savefig(os.path.join(PAP, 'fig_akw.pdf')); fig.savefig(os.path.join(PAP, 'fig_akw.png'), dpi=200)
    plt.close(fig); print('saved fig_akw ->', PAP)


if __name__ == '__main__':
    main()
