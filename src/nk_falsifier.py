# -*- coding: utf-8 -*-
"""The momentum-distribution falsifier n_k: a DIAGONAL, post-selection-protected, DISCRIMINATING
on-device channel -- the sweet spot the off-diagonal m1(k) misses.

n_k = int_{omega<0} A(k,omega) domega = <c^dag_{k,sigma} c_{k,sigma}> is the occupied momentum weight.
Unlike the total-weight sum rule m0=<{c,c^dag}>=1 (pinned -> blind) it is UNPINNED and carries a Fermi
step (strongly k-dependent -> discriminating). Unlike the first moment m1(k) (needs the hopping term ->
off-diagonal -> device-infeasible, 37-115% depolarizing bias) it is a NUMBER operator in the momentum
basis: measured by a number-conserving fermionic Fourier transform (Givens network) + occupation readout,
DIAGONAL and protected by N,S_z post-selection -> the leading depolarizing bias is removed.

Classical gate: n_k(k) has a large Fermi step (structure ~0.8-0.9) -> discriminating.
Device feasibility (in-sector depolarizing + post-selection model, validated in the off-diagonal study):
bias(k)=p_eff*|filling-n_k| with p_eff=1-(1-eps)^{n_CZ^FT}; the ~60-CZ FT contributes ~9% of the
Fermi step (INCREMENTAL cost only -- EXCLUDES the deep state-prep depolarization, which n_k inherits
like any device observable). n_k does NOT beat the off-diagonal m1 floor at present depth; it removes the
off-diagonal-specific penalties (spreading, post-selection loss), relocating feasibility to the generic
prep-depth cost. SIM/forecast: n_k has not been run on ibm_fez.
"""
import os, sys, json, numpy as np, scipy.sparse as sp
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE); import spectral_lanczos as sl
RES=os.path.normpath(os.path.join(HERE,'..','06_results'))

def cdc(L,n,p,q):
    S,idx=sl.strings(L,n); D=len(S); rows=[];cols=[];val=[]
    for a,m in enumerate(S):
        if (m>>q)&1 and (p==q or not (m>>p)&1):
            m2=(m & ~(1<<q))|(1<<p)
            if p==q: s=1.0
            else:
                lo,hi=min(p,q),max(p,q); mask=m&(((1<<hi)-1)^((1<<(lo+1))-1)); s=-1.0 if(bin(mask).count('1')&1)else 1.0
            rows.append(idx[m2]);cols.append(a);val.append(s)
    return sp.csr_matrix((val,(rows,cols)),shape=(D,D))

def nk_true(L,U,fn):
    nup=nd=fn; Tu,Su,_=sl.hop(L,nup); Td,Sd,_=sl.hop(L,nd); Du,Dd=len(Su),len(Sd)
    up=sl.occ_matrix(Su,L); dn=sl.occ_matrix(Sd,L)
    Hs=(sp.kron(Tu,sp.identity(Dd))+sp.kron(sp.identity(Du),Td)+sp.diags(U*(up@dn.T).ravel())).tocsr()
    E0,psi=sl._sub_gs(Hs); psi=np.real(psi)
    rho1=np.array([[float(psi@(sp.kron(cdc(L,nup,p,q),sp.identity(Dd))@psi)) for q in range(L)] for p in range(L)])
    ks=np.array([2*np.pi*mm/L for mm in range(L)])
    nk=np.array([np.real(sum(np.exp(-1j*k*(p-q))*rho1[p,q] for p in range(L) for q in range(L)))/L for k in ks])
    return ks,nk,fn/L

def feasibility(L,nk,fill,eps=2.5e-3):
    nGivens=L*(L-1)//2; nCZ=nGivens*2*2                # Clements L(L-1)/2 Givens/spin, ~2 CZ each, 2 spins
    peff=1-(1-eps)**nCZ
    bias=peff*np.abs(fill-nk)
    rungs={}
    for lam in [0.10,0.20,0.35,0.50]:
        dn=lam*np.abs(nk-fill); ratio=float(bias.max()/dn.max())
        rungs[f'{lam:.2f}']={'max_Delta_n':float(dn.max()),'bias_over_signal':ratio,
                             'verdict':'feasible' if ratio<0.5 else 'marginal' if ratio<1 else 'dead'}
    return nCZ,float(peff),float(bias.max()),rungs

out={'_provenance':{'script':'nk_falsifier.py','sim_only':True,
     'model':'in-sector depolarizing + N,Sz post-selection (validated in off-diagonal study); FT CZ from Clements decomposition',
     'note':'device-feasibility FORECAST; n_k not yet run on ibm_fez'}}
for L,U,fn in [(6,4.0,2),(6,8.0,2),(6,8.0,3)]:
    ks,nk,fill=nk_true(L,U,fn)
    nCZ,peff,biasmax,rungs=feasibility(L,nk,fill)
    out[f'L{L}_U{U}_n{fn}']={'filling':fill,'ks_over_pi':[float(k/np.pi) for k in ks],'n_k':[float(v) for v in nk],
        'structure_maxmin':float(nk.max()-nk.min()),'FT_CZ':nCZ,'p_eff':peff,'bias_max':biasmax,'feasibility':rungs}
    print(f"L={L} U/t={U} n={fn}: n_k structure={nk.max()-nk.min():.3f}  FT~{nCZ}CZ p_eff={peff:.3f} bias_max={biasmax:.3f}")
    for lam,r in rungs.items(): print(f"    lam={lam}: Dn={r['max_Delta_n']:.3f} bias/sig={r['bias_over_signal']:.2f} {r['verdict']}")
json.dump(out,open(os.path.join(RES,'2026-09-01_nk_falsifier.json'),'w'),indent=2)
print("\nwrote 06_results/2026-09-01_nk_falsifier.json")
