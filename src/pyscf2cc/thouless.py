import os
import numpy as np
from pyscf.ao2mo.incore import general

def thouless_orb(orb, f_T1 = 'T1.dat'):
    '''
    get Thouless orbital: |\Tilde{i}>=|i>+t(a,i)|a>
    t(a,i) is CC T1 amplitudes
    the Thouless orbital form a non-orthogonal determinant |\phi>=e^T1 |0>

    orb: orbital object
    f_T1: binary file to store (T1AA,T1BB), need store as (i,a), and without symmerty

    return:
        Thouless orbitals after normalization.
    '''
    mf = orb.mf
    mol = orb.mol
    ref = orb.ref
    if ref == 'RHF':
        coa = orb.mo_coeff[:, 0:orb.n_occ]
        cva = orb.mo_coeff[:, orb.n_occ:]
        noa = orb.n_occ
        nva = orb.n_vir
    elif ref == 'UHF':
        coa = orb.mo_coeff[0][:, 0:orb.n_occ[0]]
        cva = orb.mo_coeff[0][:, orb.n_occ[0]:]
        cob = orb.mo_coeff[1][:, 0:orb.n_occ[1]]
        cvb = orb.mo_coeff[1][:, orb.n_occ[1]:]
        noa = orb.n_occ[0]
        nob = orb.n_occ[1]
        nva = orb.n_vir[0]
        nvb = orb.n_vir[1]

    #------ get T1
    number = 10
    t1AA = np.fromfile(file=orb.tmpdir+f_T1, dtype=np.float64, count=noa*nva, offset=0) # t(i,a)
    t1AA = t1AA.reshape(noa,nva)
    idx = np.argsort(abs(t1AA).ravel())[-number:]
    row, col = np.unravel_index(idx, t1AA.shape)
    print('T1AA:')
    for i in range(0,number):
        j = number-1-i
        print(row[j],'--->', col[j], t1AA[row[j],col[j]])

    if ref == 'UHF':
        t1BB = np.fromfile(file=orb.tmpdir+f_T1, dtype=np.float64, count=nob*nvb, offset=8*noa*nva) # every float64 with 8 byte
        t1BB = t1BB.reshape(nob,nvb)
        idx = np.argsort(abs(t1BB).ravel())[-number:]
        row, col = np.unravel_index(idx, t1BB.shape)
        print('T1BB:')
        for i in range(0,number):
            j = number-1-i
            print(row[j],'--->', col[j], t1BB[row[j],col[j]])

    #------ get Thouless orbitals
    print('Alpha Thouless orbital')
    ovlp = mf.get_ovlp(mol)
    # \tilde{c}(u,i)=c(u,i)+t(i,a)c(u,a)
    Thouless_orb_a = np.einsum('ua,ia->ui', cva, t1AA)
    Thouless_orb_a += coa
    Thouless_orb_a = _normalize(Thouless_orb_a, ovlp)
    if ref == 'RHF':
        return Thouless_orb_a

    print('Beta Thouless orbital')
    Thouless_orb_b = np.einsum('ua,ia->ui', cvb, t1BB)
    Thouless_orb_b += cob
    Thouless_orb_b = _normalize(Thouless_orb_b, ovlp)
    return [Thouless_orb_a,Thouless_orb_b]

def _normalize(coeff, ovlp):
    '''
    Normalize the Thouless orbital and determinant
    '''
    norm, ovlp_mo = det_ovlp(coeff, coeff, ovlp)
    norm = np.sqrt(norm)
    print('Thouless orbital ovlp matrix:\n', ovlp_mo)
    print('Norm of Thouless determinant:', norm)

    # normalization Thouless orbital
    nmo = ovlp_mo.shape[0]
    for i in range(0,nmo):
        coeff[:,i] = coeff[:,i]/np.sqrt(ovlp_mo[i,i])
    
    # normaliza Thouless determinant
    # for determinant state \psi:
    # <\phi|\phi> = |S|
    # where S_ij = <i|j>, i, j is the orbital form \phi
    norm, ovlp_mo = det_ovlp(coeff, coeff, ovlp)
    norm = np.sqrt(norm)
    fact = norm**(1/nmo)
    coeff = coeff/fact

    # <p|q>=c^*(u,p) c(v,q) <u|v>
    norm, ovlp_mo = det_ovlp(coeff, coeff, ovlp)
    norm = np.sqrt(norm)
    print('Thouless orbital ovlp matrix after normilization:\n', ovlp_mo)
    print('Norm of Thouless determinant after normilization:', norm)

    return coeff

def det_ovlp(coeff1, coeff2, ovlp):
    '''
    Calculate overlop of two determinant overlop:\n
    <\Tilde{\psi}|\phi> = |S|\n
    where S_ij = <\Tilde{i}|j>

    coeff1: orbital of \Tilde{\psi}
    coeff2: orbital of \phi
    ovlp: overlop matrix of AO

    return:\n
        <\Tilde{\psi}|\phi>, MO overlop matrix
    '''
    # <p|q>=c^*(u,p) c(v,q) <u|v>
    ovlp_mo = np.einsum('up,vq,uv->pq', coeff1.conj(), coeff2, ovlp)
    ovlp_det = np.linalg.det(ovlp_mo)
    return ovlp_det, ovlp_mo

def thouless_T1(current_coeff, current_occ, target_coeff, target_occ, ovlp, number=10):
    '''
    Calculate Thouless T1 form current orbital |p> to target orbital |\Tilde{m}>:\n
    t(a,i) = C(a,m) U(m,i)\n
    C(p,m) = <p|\Tilde{m}>, devide into cooupied and virtual block\n
    [U(m,i)] = [C(i,m)]^(-1)\n

    ovlp: AO overlop matrix

    return:
        t(i,a)
    '''
    current = [None]*2
    idx = [i for i,value in enumerate(current_occ) if value==1]
    current[0] = current_coeff[:,idx]

    idx = [i for i,value in enumerate(current_occ) if value==0]
    current[1] = current_coeff[:,idx]

    idx = [i for i,value in enumerate(target_occ) if value==1]
    target = target_coeff[:,idx]
    
    # calculate C(p,m) = c^*(u,p)<u|v>\tilde{c}(v,m)
    Cim = np.einsum('ui,uv,vm->im',current[0].conj(), ovlp, target)
    Cam = np.einsum('ua,uv,vm->am',current[1].conj(), ovlp, target)

    # calculate inverse of Cim
    Umi = np.linalg.inv(Cim)
    
    # calculate t(i,a)
    tia = np.einsum('mi,am->ia', Umi, Cam)

    number = 10
    idx = np.argsort(abs(tia).ravel())[-number:]
    row, col = np.unravel_index(idx, tia.shape)
    print('Thouless T1:')
    for i in range(0,number):
        j = number-1-i
        print(row[j],'--->', col[j], tia[row[j],col[j]])

    return tia

def thouless_T1_det(current_coeff, current_occ, target, ovlp):
    '''
    Calculate Thouless T1 form current orbital |p> to target determinant |\Tilde{m}>:\n
    t(a,i) = C(a,m) U(m,i)\n
    C(p,m) = <p|\Tilde{m}>, devide into cooupied and virtual block\n
    [U(m,i)] = [C(i,m)]^(-1)\n

    target: target determinant |\Tilde{m}>
    ovlp: AO overlop matrix

    return:
        t(i,a)
    '''
    current = [None]*2
    idx = [i for i,value in enumerate(current_occ) if value==1]
    current[0] = current_coeff[:,idx]

    idx = [i for i,value in enumerate(current_occ) if value==0]
    current[1] = current_coeff[:,idx]

    # calculate C(p,m) = c^*(u,p)<u|v>\tilde{c}(v,m)
    Cim = np.einsum('ui,uv,vm->im',current[0].conj(), ovlp, target)
    Cam = np.einsum('ua,uv,vm->am',current[1].conj(), ovlp, target)

    # calculate inverse of Cim
    Umi = np.linalg.inv(Cim)
    
    # calculate t(i,a)
    tia = np.einsum('mi,am->ia', Umi, Cam)

    number = 10
    idx = np.argsort(abs(tia).ravel())[-number:]
    row, col = np.unravel_index(idx, tia.shape)
    print('Thouless T1:')
    for i in range(0,number):
        j = number-1-i
        print(row[j],'--->', col[j], tia[row[j],col[j]])

    return tia

def thouless_orb_incore(coeff, T1, ovlp):
    '''
    get Thouless orbital: |\Tilde{i}>=|i>+t(i,a)|a>
    t(i,a) is CC T1 amplitudes
    the Thouless orbital form a non-orthogonal determinant |\phi>=e^T1 |0>

    coeff: list, [|i>,|a>] 
    T1: t(i,a)
    ovlp: AO overlop matrix

    return:
        Thouless orbitals after normalization.
    '''
    # \tilde{c}(u,i)=c(u,i)+t(i,a)c(u,a)
    Thouless_orb = np.einsum('ua,ia->ui', coeff[1], T1)
    Thouless_orb += coeff[0]
    Thouless_orb = _normalize(Thouless_orb, ovlp)

    return Thouless_orb

def save_t1(t1, file='T1.dat', path='./tmp/'):
    '''
    Save flat T1 in file:
    '''
    if (not isinstance(t1, list)) and np.iscomplexobj(t1):
        # (t1.real, t1.imag)
        print('WARN: GHF complex t1')
        data = np.concatenate([t1.ravel().real, t1.ravel().imag])
    else:
        # (t1AA,t1BB)
        data = np.concatenate([t1[0].ravel(), t1[1].ravel()])
    
    os.makedirs(path, exist_ok=True)
    data.tofile(path+file)

def select_orbital(coeff1, iorb1, coeff2, ovlp):
    '''
    For an virtual orbital in set1, select an virtual orbital in set 2,\n
    to make the <det1|det2> have max overlop

    coeff1: list, [occ1,vir1]
    iorb1: list, idx of virtual orbital we need in vir1
    coeff2: list, [occ2,vir2]
    ovlp: AO overlop matrix
    
    return:
        iorb2: list, idx of virtual orbital in vir2 correspond to iorb1
    '''
    occ1 = coeff1[0]
    vir1 = coeff1[1]
    occ2 = coeff2[0]
    vir2 = coeff2[1]
    target1 = vir1[:,iorb1] # virtual orbital we need in vir1

    ovlp_mo = np.einsum('up,vq,uv->pq', target1.conj(), vir2, ovlp)
    iorb2 = []
    # cycle for orbital in vir1
    for i in range(0,len(iorb1)):
        print('--- i:', iorb1[i])
        # First, I use the orbital overlop to select the most important vir2 for targer1
        ovlp_i = ovlp_mo[i,:]
#        idx = np.where(np.abs(ovlp_i) > 0.1)[0] # index of some max overlop orbital
        idx = np.argpartition(np.abs(ovlp_i), -3)[-3:]

        orb1 = np.hstack((occ1, target1[:,i].reshape(-1,1))) # det1
        ovlp_det = []
        for j in range(0,len(idx)):
            # Second, I calculate <det1|det2> to select the best orbital
            orb2 = np.hstack((occ2, vir2[:,idx[j]].reshape(-1,1))) # det2
            ovlp_det.append(det_ovlp(orb1, orb2, ovlp)[0])

            print('j:', idx[j], ovlp_i[idx[j]], ovlp_det[j])

        max_j = np.argmax(np.abs(ovlp_det))
        iorb2.append(idx[max_j])

    return iorb2

def select_orbital_det(coeff, coeff_target, ovlp):
    '''
    For an determinant det2, select an virtual orbital in set 1,\n
    to make the <det1|det2> have max overlop

    coeff: list, [occ1,vir1]
    coeff_target: orbital formed det2
    ovlp: AO overlop matrix
    
    return:
        max_i: idx of virtual orbital in vir1
    '''
    occ1 = coeff[0]
    vir1 = coeff[1]

    # cycle for orbital in vir1
    ovlp_det = []
    for i in range(0,vir1.shape[1]):
        orb1 = np.hstack((occ1, vir1[:,i].reshape(-1,1))) # det1
        # I calculate <det1|det2> to select the best orbital
        ovlp_det.append(det_ovlp(orb1, coeff_target, ovlp)[0])

        print('i:', i, ovlp_det[i])

    max_i = np.argmax(np.abs(ovlp_det))

    return max_i

def spin_coupling(orb, coeff, h1e_ao=None, v2e_ao=None, ovlp_ao=None):
    '''
    Calculate <tilde_Phi|H|Phi>, tilde_Phi is spin flip of Phi

    orb: orbital object
    coeff: Phi, [alpha, beta]
    '''
    mf = orb.mf
    mol = orb.mol
    if h1e_ao is None:
        h1e_ao = mf.get_hcore()
    if v2e_ao is None:
        v2e_ao = mol.intor('int2e', aosym='s8')
    if ovlp_ao is None:
        ovlp_ao = mf.get_ovlp(mol)

    # get biorthogonal MO
    l_alpha, r_alpha, D_alpha = biorthogonal(coeff[1], coeff[0], ovlp_ao)
    l_beta = r_alpha.copy()
    r_beta = l_alpha.copy()
    D_beta = D_alpha.copy()
    dim = len(D_alpha)
    print('Singular value:\n',D_alpha,'\n',D_beta)

    # get h1e_ij = <tilde_i|h|j>
    h1e_aa = np.einsum('uv,ui,vj->ij', h1e_ao, l_alpha, r_alpha)
    h1e_bb = np.einsum('uv,ui,vj->ij', h1e_ao, l_beta, r_beta)

    # get v2e_ijkl = (tilde_i,j|V|tulde_k,l)
    v2e_aaaa = general(v2e_ao, [l_alpha,r_alpha,l_alpha,r_alpha], compact=False)
    v2e_aaaa = v2e_aaaa.reshape(dim,dim,dim,dim)
    v2e_bbbb = general(v2e_ao, [l_beta,r_beta,l_beta,r_beta], compact=False)
    v2e_bbbb = v2e_bbbb.reshape(dim,dim,dim,dim)
    v2e_abab = general(v2e_ao, [l_alpha,r_alpha,l_beta,r_beta], compact=False)
    v2e_abab = v2e_abab.reshape(dim,dim,dim,dim)

    # overlop of alpha and beta part
    S_alpha = np.prod(D_alpha)
    S_beta = np.prod(D_beta)

    # H1e = \sum_i^n{D_1...h_ii...D_n}
    H1e = 0.0
    def term1(D,h1e,i):
        tmp = D[i]
        D[i] = h1e[i,i]
        value = np.prod(D)
        D[i] = tmp
        return value
    # alpha h
    for i in range(0,dim):
        H1e += term1(D_alpha,h1e_aa,i)*S_beta
    # beta h
    for i in range(0,dim):
        H1e += term1(D_beta,h1e_bb,i)*S_alpha
    print('H1e:',H1e)

    # V2e = \sum_(i<j){D_1...(V_iijj-V_ijji)...D_n}
    V2e = 0.0
    def term2(D,v2e,i,j):
        tmp1 = D[i]
        tmp2 = D[j]
        D[i] = 1.0
        D[j] = v2e[i,i,j,j]-v2e[i,j,j,i]
        value = np.prod(D)
        D[i] = tmp1
        D[j] = tmp2
        return value
    # AAAA
    for i in range(0,dim):
        for j in range(i,dim):
            V2e += term2(D_alpha,v2e_aaaa,i,j)*S_beta
    print('V2e:', V2e)
    # BBBB
    for i in range(0,dim):
        for j in range(i,dim):
            V2e += term2(D_beta,v2e_bbbb,i,j)*S_alpha
    print('V2e:', V2e)
    # ABAB
    def term3(D_a,D_b,v2e,i,j):
        tmp1 = D_a[i]
        tmp2 = D_b[j]
        D_a[i] = 1.0
        D_b[j] = v2e[i,i,j,j]
        value = np.prod(D_a)*np.prod(D_b)
        D_a[i] = tmp1
        D_b[j] = tmp2
        return value
    for i in range(0,dim):
        for j in range(0,dim):
            V2e += term3(D_alpha,D_beta,v2e_abab,i,j)
    print('V2e:', V2e)

    print('H1e+V2e:', H1e+V2e)
    return H1e+V2e

def biorthogonal(coeff_l, coeff_r, ovlp_ao):
    '''
    Biorthogonal left and right MO use SVD
    new MO satisfy <L_i|R_j> = delta(i,j)*D_i

    retrun: new_l, new_r, singular_value
    '''

    # S_ij = <L_i|R_j>
    S = np.einsum('uv,up,vq->pq', ovlp_ao, coeff_l, coeff_r)
    # S = U*D*V^T, U and V are orthogonal matrix, D is diagonal matrix
    U, D, VT = np.linalg.svd(S)
    V = VT.transpose(1,0)
    # calculate new basis <tilde_i|=<j|U_ji and |tilde_i>=VT_ij|j>
    # which satisfy <tilde_i|tilde_j>=D_ij
    new_l = np.einsum('uj,ji->ui', coeff_l, U)
    new_r = np.einsum('uj,ji->ui', coeff_r, V)

    return new_l, new_r, D