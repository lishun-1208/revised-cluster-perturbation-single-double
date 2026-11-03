from pathlib import Path
import numpy as np
import pyscf2cc.symm_info as si
from pyscf import scf, gto, ao2mo, lib
import pyscf2cc.write_int as wi
import pyscf2cc.interface as inter
from time import process_time, time
import h5py

def py2cc(mf_in, mol, ref, fz=0, mo_occ=None, canonical=False, test=False, to_cc=True):
    '''
    mf_in: scf or casscf object, only support RCASSCF
    mol: gto.M()
    ref: str, UHF, RHF, ROHF, CASSCF
    fz: int, the number of frozen orbital
    mo_occ: 2D array, orbital occupation
    canonical: bool, whether canonicalize orbital
    test: bool, whether run in test module
    to_cc: bool, whether save file for later CC calculation

    return:
        orb: orbital object
    '''

    # convert to scf object which can use for symmerty class
    if ref == 'UHF':
        mf = scf_copy(mf_in)
        if mo_occ is not None:
            mf.mo_occ = mo_occ

    elif ref == 'RHF':
        mf = scf_copy(mf_in)
        if mo_occ is not None:
            if not np.array_equal(mo_occ[0], mo_occ[1]):
                print('convert RHF to UHF')
                mf = scf.addons.convert_to_uhf(mf)
                mf.mo_occ = mo_occ
                ref = 'UHF'
            else:
                mf.mo_occ = np.array(mo_occ[0])*2

    elif ref == 'ROHF':
        mf = scf_copy(mf_in)
        if mo_occ is None:
            print('convert ROHF to UHF')
            mf = scf.addons.convert_to_uhf(mf)
            ref = 'UHF'
        else:
            if not np.array_equal(mo_occ[0], mo_occ[1]):
                print('convert ROHF to UHF')
                mf = scf.addons.convert_to_uhf(mf)
                mf.mo_occ = mo_occ
                ref = 'UHF'
            else:
                print('convert ROHF to RHF')
                print('WARN: will create a new SCF object')
                ref = 'RHF'
                mf2 = scf.RHF(mol)
                mf2.mo_coeff = mf.mo_coeff
                mf2.mo_energy = mf.mo_energy
                mf2.mo_occ = np.array(mo_occ[0])*2
                mf = mf2

    elif ref == 'CASSCF':
        if mo_occ is None:
            print('ERROR: to use CASSCF as reference, must give a mo_occ')
            exit()

        if np.array_equal(mo_occ[0], mo_occ[1]):
            print('convert CASSCF to RHF')
            print('WARN: will create a new SCF object')
            ref = 'RHF'
            mf2 = scf.RHF(mol)
            mf2.mo_coeff = mf.mo_coeff
            mf2.mo_energy = mf.mo_energy
            mf2.mo_occ = np.array(mo_occ[0])*2
        else:
            print('convert CASSCF to UHF')
            print('WARN: will create a new SCF object')
            ref = 'UHF'
            mf2 = scf.UHF(mol)
            mf2.mo_coeff = [None]*2
            mf2.mo_coeff[0] = mf.mo_coeff.copy()
            mf2.mo_coeff[1] = mf.mo_coeff.copy()
            mf2.mo_energy = [None]*2
            mf2.mo_energy[0] = mf.mo_energy.copy()
            mf2.mo_energy[1] = mf.mo_energy.copy()
            mf2.mo_occ = np.array(mo_occ)
        mf = mf2

    tmpdir = _read_tmpdir()
    if not Path(tmpdir).exists():
        Path(tmpdir).mkdir(parents=True, exist_ok=True)

    if canonical:
        mf = si.canonical(mf)

    print('mf class name:',mf.__class__.__name__)
    orb = si.orbital(mf, mol, fz, ref, test, tmpdir)
    # get orbital information
    orb.get_sym_info()
    orb.resort_orb()

    if to_cc:
        print('Save imtermidate file for CC calculation.')
        write_interface(orb)
        get_2e_int_outcore(orb, tmpdir)
        get_fock(orb, file_path=tmpdir)

    return orb

def get_2e_int(orb, integral='int2e', file_path=None, file_head=' ', sym_A=0):
    '''
    orb: orbital object
    integral: str, the name of integral need get, see pyscf/gto/moleintor.py
    sym_A: symmerty of operator A
    file_head: str with len=1

    Get 2e integral A^pq_rs = <pq|A|rs>-<pq|A|sr> and write in hard disk in our own order.
    PySCF store 2e-integral as (pq|rs)=<pr|qs>
    '''

    if file_path is None:
        file_path = orb.tmpdir

    mol = orb.mol
    ref = orb.ref

    # get orbital coefficient
    if ref == 'RHF':
        coeff_occ_a = orb.mo_coeff[:, 0:orb.n_occ]
        coeff_vir_a = orb.mo_coeff[:, orb.n_occ:]
        coeff_occ_b = coeff_occ_a
        coeff_vir_b = coeff_vir_a
    elif ref == 'UHF':
        coeff_occ_a = orb.mo_coeff[0][:, 0:orb.n_occ[0]]
        coeff_vir_a = orb.mo_coeff[0][:, orb.n_occ[0]:]
        coeff_occ_b = orb.mo_coeff[1][:, 0:orb.n_occ[1]]
        coeff_vir_b = orb.mo_coeff[1][:, orb.n_occ[1]:]
    print('shape of coeff_occ_a', coeff_occ_a.shape)
    print('shape of coeff_vir_a', coeff_vir_a.shape)
    print('shape of coeff_occ_b', coeff_occ_b.shape)
    print('shape of coeff_vir_b', coeff_vir_b.shape)

    # get AO integral (pq|rs)
    ao_int = mol.intor(integral, aosym='s1')
    print('2e ao int shape',ao_int.shape)
    if ao_int.size != orb.norb_no_fz**4:
        print('ERROR: size of ao integral error', ao_int.size, orb.norb_no_fz**4)
    # resort: (pq|rs)-><pr|qs>
    ao_int = ao_int.transpose(0,2,1,3)
#    ao_int = mol.intor(integral, aosym='s8')

    # initialize mod
    wi.init(orb.dirprd+1, orb.nirrep, file_path)

    # get integral for different spin and orbital type
    _get_2e_spin_int(orb, ao_int, coeff_occ_a, coeff_occ_b, coeff_vir_a, coeff_vir_b, \
                     sym_A, file_head)

def _get_2e_spin_int(orb, ao_int, c_o_a, c_o_b, c_v_a, c_v_b, \
                     sym_A, file_head):
    '''
    c_o/v_a/b: mo coefficient for occupied and virtual spin alpha and beta orbital
    '''
    ref = orb.ref
    popa = orb.popa
    popb = orb.popb
    vira = orb.vira
    virb = orb.virb

    nCore = 100000000
    icore = np.empty(nCore, dtype=np.float64) # 800MB

    if file_head == ' ':
        default_name = True
    else:
        default_name = False

    #------ get abcd integral
    t1 = process_time()
    wt1 = time()
    # ABAB
    if default_name:
        file = 'ABCDAB__'
    else:
        file = file_head + '2e1ABAB'
    _2e_mo_int(orb, ao_int, c_v_a, c_v_b, c_v_a, c_v_b, \
               vira, virb, vira, virb, sym_A, False, False, icore, file, 'vvvv')
    if ref == 'UHF':
        # AAAA
        if default_name:
            file = 'ABCDAA__'
        else:
            file = file_head + '2e1AAAA'
        _2e_mo_int(orb, ao_int, c_v_a, c_v_a, c_v_a, c_v_a, \
                   vira, vira, vira, vira, sym_A, True, True, icore, file, 'vvvv')
        # BBBB
        if default_name:
            file = 'ABCDBB__'
        else:
            file = file_head + '2e1BBBB'
        _2e_mo_int(orb, ao_int, c_v_b, c_v_b, c_v_b, c_v_b, \
                   virb, virb, virb, virb, sym_A, True, True, icore, file, 'vvvv')
    t2 = process_time()
    wt2 = time()
    print('cpu time for get VVVV integral:', t2-t1)
    print('wall time for get VVVV integral:', wt2-wt1)

    #------ get abij integral
    t1 = process_time()
    wt1 = time()
    # ABAB
    if default_name:
        file = 'ABIJAB__'
    else:
        file = file_head + '2e2ABAB'
    _2e_mo_int(orb, ao_int, c_v_a, c_v_b, c_o_a, c_o_b, \
               vira, virb, popa, popb, sym_A, False, False, icore, file, 'vvoo')
    if ref == 'UHF':
        # AAAA
        if default_name:
            file = 'ABIJAA__'
        else:
            file = file_head + '2e2AAAA'
        _2e_mo_int(orb, ao_int, c_v_a, c_v_a, c_o_a, c_o_a, \
                   vira, vira, popa, popa, sym_A, True, True, icore, file, 'vvoo')
        # BBBB
        if default_name:
            file = 'ABIJBB__'
        else:
            file = file_head + '2e2BBBB'
        _2e_mo_int(orb, ao_int, c_v_b, c_v_b, c_o_b, c_o_b, \
                   virb, virb, popb, popb, sym_A, True, True, icore, file, 'vvoo')
    t2 = process_time()
    wt2 = time()
    print('cpu time for get VVOO integral:', t2-t1)
    print('wall time for get VVOO integral:', wt2-wt1)

    #------ get ijkl integral
    t1 = process_time()
    wt1 = time()
    # ABAB
    if default_name:
        file = 'IJKLAB__'
    else:
        file = file_head + '2e3ABAB'
    _2e_mo_int(orb, ao_int, c_o_a, c_o_b, c_o_a, c_o_b, \
               popa, popb, popa, popb, sym_A, False, False, icore, file, 'oooo')
    if ref == 'UHF':
        # AAAA
        if default_name:
            file = 'IJKLAA__'
        else:
            file = file_head + '2e3AAAA'
        _2e_mo_int(orb, ao_int, c_o_a, c_o_a, c_o_a, c_o_a, \
                   popa, popa, popa, popa, sym_A, True, True, icore, file, 'oooo')
        # BBBB
        if default_name:
            file = 'IJKLBB__'
        else:
            file = file_head + '2e3BBBB'
        _2e_mo_int(orb, ao_int, c_o_b, c_o_b, c_o_b, c_o_b, \
                   popb, popb, popb, popb, sym_A, True, True, icore, file, 'oooo')
    t2 = process_time()
    wt2 = time()
    print('cpu time for get OOOO integral:', t2-t1)
    print('wall time for get OOOO integral:', wt2-wt1)

    #------ get abci integral
    t1 = process_time()
    wt1 = time()
    # ABAB
    if default_name:
        file = 'ABCIAB__'
    else:
        file = file_head + '2e4ABAB'
    _2e_mo_int(orb, ao_int, c_v_a, c_v_b, c_v_a, c_o_b, \
               vira, virb, vira, popb, sym_A, False, False, icore, file, 'vvvo')
    if ref == 'UHF':
        # AAAA
        if default_name:
            file = 'ABCIAA__'
        else:
            file = file_head + '2e4AAAA'
        _2e_mo_int(orb, ao_int, c_v_a, c_v_a, c_v_a, c_o_a, \
                   vira, vira, vira, popa, sym_A, True, False, icore, file, 'vvvo')
        # BBBB
        if default_name:
            file = 'ABCIBB__'
        else:
            file = file_head + '2e4BBBB'
        _2e_mo_int(orb, ao_int, c_v_b, c_v_b, c_v_b, c_o_b, \
                   virb, virb, virb, popb, sym_A, True, False, icore, file, 'vvvo')
        # BABA
        if default_name:
            file = 'ABCIBA__'
        else:
            file = file_head + '2e4BABA'
        _2e_mo_int(orb, ao_int, c_v_b, c_v_a, c_v_b, c_o_a, \
                   virb, vira, virb, popa, sym_A, False, False, icore, file, 'vvvo')
    t2 = process_time()
    wt2 = time()
    print('cpu time for get VVVO integral:', t2-t1)
    print('wall time for get VVVO integral:', wt2-wt1)

    #------ get ijka integral
    t1 = process_time()
    wt1 = time()
    # ABAB
    if default_name:
        file = 'IJKAAB__'
    else:
        file = file_head + '2e5ABAB'
    _2e_mo_int(orb, ao_int, c_o_a, c_o_b, c_o_a, c_v_b, \
               popa, popb, popa, virb, sym_A, False, False, icore, file, 'ooov')
    if ref == 'UHF':
        # AAAA
        if default_name:
            file = 'IJKAAA__'
        else:
            file = file_head + '2e5AAAA'
        _2e_mo_int(orb, ao_int, c_o_a, c_o_a, c_o_a, c_v_a, \
                   popa, popa, popa, vira, sym_A, True, False, icore, file, 'ooov')
        # BBBB
        if default_name:
            file = 'IJKABB__'
        else:
            file = file_head + '2e5BBBB'
        _2e_mo_int(orb, ao_int, c_o_b, c_o_b, c_o_b, c_v_b, \
                   popb, popb, popb, virb,  sym_A, True, False, icore, file, 'ooov')
        # BABA
        if default_name:
            file = 'IJKABA__'
        else:
            file = file_head + '2e5BABA'
        _2e_mo_int(orb, ao_int, c_o_b, c_o_a, c_o_b, c_v_a, \
                   popb, popa, popb, vira, sym_A, False, False, icore, file, 'ooov')
    t2 = process_time()
    wt2 = time()
    print('cpu time for get OOOV integral:', t2-t1)
    print('wall time for get OOOV integral:', wt2-wt1)

    #------ get iabj integral
    t1 = process_time()
    wt1 = time()
    # AAAA
    if default_name:
        file = 'BIAJAAAA'
    else:
        file = file_head + '2e6AAAA'
    _2e_mo_int_iabj(orb, ao_int, c_o_a, c_v_a, c_v_a, c_o_a, \
                    popa, vira, vira, popa, sym_A, icore, file, 1)
    # ABAB
    if default_name:
        file = 'BIAJAABB'
    else:
        file = file_head + '2e6AABB'
    _2e_mo_int_iabj(orb, ao_int, c_o_a, c_v_b, c_v_a, c_o_b, \
                    popa, virb, vira, popb, sym_A, icore, file, 2)
    # BAAB
    if default_name:
        file = 'BIAJABAB'
    else:
        file = file_head + '2e6ABAB'
    _2e_mo_int_iabj(orb, ao_int, c_o_b, c_v_a, c_v_a, c_o_b, \
                    popb, vira, vira, popb, sym_A, icore, file, 3)
    if ref == 'UHF':
        # BBBB
        if default_name:
            file = 'BIAJBBBB'
        else:
            file = file_head + '2e6BBBB'
        _2e_mo_int_iabj(orb, ao_int, c_o_b, c_v_b, c_v_b, c_o_b, \
                        popb, virb, virb, popb, sym_A, icore, file, 1)
        # ABBA
        if default_name:
            file = 'BIAJBABA'
        else:
            file = file_head + '2e6BABA'
        _2e_mo_int_iabj(orb, ao_int, c_o_a, c_v_b, c_v_b, c_o_a, \
                        popa, virb, virb, popa, sym_A, icore, file, 3)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get OVVO integral:', t2-t1)
    print('wall time for get OVVO integral:', wt2-wt1)

def _2e_mo_int(orb, ao_int, c_p, c_q, c_r, c_s, \
               pop_p, pop_q, pop_r, pop_s, sym_I, \
               anti_pq, anti_rs, icore, file, type):
    '''
    anti_pq: whether anti symmetrize pq
    anti_rs: whether anti symmetrize rs
    c_p,q,r,s: coefficient
    type: 'vvvv', 'vvoo', 'oooo', 'vvvo', 'ooov'
    
    Calculate MO interal I(pq,rs) for abcd, abij, ijkl, abci, ijka
    '''

    # create file
    if anti_pq:
        ioppq = 0
    else:
        ioppq = 1
    if anti_rs:
        ioprs = 0
    else:
        ioprs = 1
    wi.creat_file(sym_I+1, pop_p, pop_q, ioppq, pop_r, pop_s, ioprs, file, icore)

    # the scale be n^4*p+n^3*p*q+n^2*p*q*r+n*p*q*r*s
    # so its best to first contract the small p
    if type in ['oooo', 'ooov']:
        mix_int0 = lib.einsum('ijkl,ip,jq->pqkl', ao_int, c_p, c_q) # n^4*o+n^3*o^2
    
    # cycle for irrep_rs
    for irrep in range(0, orb.nirrep):
        irrep_pq = orb.dirprd[irrep, sym_I]
        n_rs = wi.calc_num_of_pq(pop_r, pop_s, irrep+1, ioprs)

        # cycle for rs
        old_s = 0
        for rs in range(0, n_rs):
            r, s = wi.get_rs(rs+1, irrep+1, pop_r, pop_s, ioprs)
            if s != old_s:
                if type in ['oooo', 'ooov']:
                    mix_int = lib.einsum('pqkl,l->pqk', mix_int0, c_s[:, s-1]) # n^2*o^2
                else:
                    mix_int = lib.einsum('ijkl,l->ijk', ao_int, c_s[:, s-1]) # n^4
                old_s = s
            _2e_mo_int_block(orb, mix_int, c_p, c_q, c_r, r, rs, \
                             pop_p, pop_q, irrep_pq, sym_I, ioppq, icore, file, type)

def _2e_mo_int_block(orb, mix_int, c_p, c_q, c_r, r, rs, \
                     pop_p, pop_q, sym_pq, sym_I, ioppq, icore, file, type):
    '''
    mix_int: mixture integral of AO and MO

    Calculate MO integral for specified r and s block
    '''

    # mo(p,q,r,s) = ao(i,j,k,l) c(ip) c(jq) c(kr) c(ls)
    # <pq|rs>=(pr|qs)
#    mo_int = ao2mo.incore.general(ao_int,[c_p, c_r[:, r-1].reshape(-1, 1),\
#                                          c_q, c_s[:, s-1].reshape(-1, 1)], compact=False)
    mo_int = lib.einsum('ijk,k->ij', mix_int, c_r[:, r-1]) # n^3
    if type not in ['oooo', 'ooov']:
        mo_int = lib.einsum('ij,ip->pj', mo_int, c_p) # n^2*(v or o)
        mo_int = lib.einsum('pj,jq->pq', mo_int, c_q) # n*(v or o)^2

    if ioppq == 0:
        # <pq||rs> = <pq|rs>-<qp|rs> = (pr|qs)-(qr|ps)
        mo_int = mo_int - mo_int.transpose(1,0)

    # change to Fortran 1D array
    mo_int = mo_int.ravel(order='F')

    # remove the zero block
    n_pq = wi.calc_num_of_pq(pop_p, pop_q, sym_pq+1, 1)
    n_p = sum(pop_p)
    n_q = sum(pop_q)
    mo_int = wi.compress_by_symmerty(n_p, n_q, n_pq, pop_p, pop_q, sym_pq+1, mo_int)

    # store in file
    wi.write_block(n_pq, mo_int, sym_pq+1, sym_I+1, pop_p, pop_q, ioppq, file, rs+1, icore)

def _2e_mo_int_iabj(orb, ao_int, c_i, c_a, c_b, c_j, \
                    pop_i, pop_a, pop_b, pop_j, sym_I, \
                    icore, file, type):
    '''
    type: 1 - AAAA, BBBB; 2 - ABAB; 3 - ABBA, BAAB
    
    Calculate MO interal iabj(bi,aj)
    '''

    # create file
    wi.creat_file(sym_I+1, pop_b, pop_i, 1, pop_a, pop_j, 1, file, icore)

    
    # cycle for irrep_aj
    for irrep in range(0, orb.nirrep):
        irrep_bi = orb.dirprd[irrep, sym_I]
        n_aj = wi.calc_num_of_pq(pop_a, pop_j, irrep+1, 1)

        # cycle for aj
        old_j = 0
        for aj in range(0, n_aj):
            a, j = wi.get_rs(aj+1, irrep+1, pop_a, pop_j, 1)
            if j != old_j:
                mix_int = lib.einsum('uvwx,x->uvw', ao_int, c_j[:, j-1])
                old_j = j
            _2e_mo_int_iabj_block(orb, mix_int, c_i, c_a, c_b, a, aj, \
                                  pop_i, pop_b, irrep_bi, sym_I, icore, file, type)

def _2e_mo_int_iabj_block(orb, mix_int, c_i, c_a, c_b, a, aj, \
                          pop_i, pop_b, sym_bi, sym_I, icore, file, type):
    '''
    type: 1 - AAAA, BBBB; 2 - ABAB; 3 - ABBA, BAAB
    
    Calculate MO interal iabj(bi,aj)
    '''

    if type == 1:
        # <ia||bj>=<ia|bj>-<ai|bj>
        mo_int = lib.einsum('uvw,v->uw', mix_int, c_a[:, a-1])-\
            lib.einsum('uvw,u->vw', mix_int, c_a[:, a-1])
    elif type == 2:
        # <iA||bJ>=<iA|bJ>
        mo_int = lib.einsum('uvw,v->uw', mix_int, c_a[:, a-1])
    elif type == 3:
        # <iA||Bj>=-<Ai|Bj>
        mo_int = -lib.einsum('uvw,u->vw', mix_int, c_a[:, a-1])
    mo_int = lib.einsum('uv,ui->iv', mo_int, c_i)
    mo_int = lib.einsum('iv,vb->ib', mo_int, c_b)
    # ib->bi
    mo_int = mo_int.transpose(1,0)

    # change to Fortran 1D array
    mo_int = mo_int.ravel(order='F')

    # remove the zero block
    n_bi = wi.calc_num_of_pq(pop_b, pop_i, sym_bi+1, 1)
    n_b = sum(pop_b)
    n_i = sum(pop_i)
    mo_int = wi.compress_by_symmerty(n_b, n_i, n_bi, pop_b, pop_i, sym_bi+1, mo_int)

    # store in file
    wi.write_block(n_bi, mo_int, sym_bi+1, sym_I+1, pop_b, pop_i, 1, file, aj+1, icore)

def get_1e_int(orb, integral, file_path=None, file_head=' ', sym_A=0, components=1):
    '''
    orb: orbital object
    integral: str, the name of integral need get, see pyscf/gto/moleintor.py
    sym_A: int or list of int, symmerty of operator A
    file_head: str or list of str with len=1
    components: components of such 1e integral

    get 1e integral <p|A|q> and write in hard disk in our own order
    '''

    if file_path is None:
        file_path = orb.tmpdir

    if file_head == ' ':
        print('ERROR: please give file_head')
        exit()

    mol = orb.mol
    ref = orb.ref

    # get orbital coefficient
    if ref == 'RHF':
        coeff_occ_a = orb.mo_coeff[:, 0:orb.n_occ]
        coeff_vir_a = orb.mo_coeff[:, orb.n_occ:]
        coeff_occ_b = coeff_occ_a
        coeff_vir_b = coeff_vir_a
    elif ref == 'UHF':
        coeff_occ_a = orb.mo_coeff[0][:, 0:orb.n_occ[0]]
        coeff_vir_a = orb.mo_coeff[0][:, orb.n_occ[0]:]
        coeff_occ_b = orb.mo_coeff[1][:, 0:orb.n_occ[1]]
        coeff_vir_b = orb.mo_coeff[1][:, orb.n_occ[1]:]
    print('shape of coeff_occ_a', coeff_occ_a.shape)
    print('shape of coeff_vir_a', coeff_vir_a.shape)
    print('shape of coeff_occ_b', coeff_occ_b.shape)
    print('shape of coeff_vir_b', coeff_vir_b.shape)

    # get AO integral <p|A|q>
    ao_int = mol.intor(integral,comp=components)
    if components == 1:
        print('1e ao int shape',ao_int.shape)
        if ao_int.size != orb.norb_no_fz**2:
            print('ERROR: size of ao integral error', ao_int.size, orb.norb_no_fz**2)
    else:
        print('1e ao int shape',ao_int[0].shape)
        if ao_int[0].size != orb.norb_no_fz**2:
            print('ERROR: size of ao integral error', ao_int[0].size, orb.norb_no_fz**2)

    # make axis transform for original to symmetry adapted
    if (orb.compGroup != 'C1') and (integral in ['int1e_ipovlp','int1e_r']):
        rota = np.array(mol._symm_axes)
        move = mol._symm_orig
        print('Rotation:\n', rota)
        print('Move:', move)
        new_ao_int = [None]*3
        for i in range(0,3):
            new_ao_int[i] = rota[i,0]*ao_int[0]+rota[i,1]*ao_int[1]+rota[i,2]*ao_int[2]
        if integral == 'int1e_r':
            print('ATTENTATION: not move coordinates origion to symmerty origion.')
        ao_int = new_ao_int

    # initialize mod
    wi.init(orb.dirprd+1, orb.nirrep, file_path)

    if components == 1:
        _get_1e_int(orb, ao_int, coeff_occ_a, coeff_occ_b, coeff_vir_a, coeff_vir_b,\
                    sym_A, file_head)
    else:
        for comp in range(0,components):
            _get_1e_int(orb, ao_int[comp], coeff_occ_a, coeff_occ_b, coeff_vir_a, coeff_vir_b,\
                        sym_A[comp], file_head[comp])

def _get_1e_int(orb, ao_int, c_o_a, c_o_b, c_v_a, c_v_b, \
                sym_A, file_head):
    '''
    c_o/v_a/b: mo coefficient for occupied and virtual spin alpha and beta orbital
    '''
    ref = orb.ref
    popa = orb.popa
    popb = orb.popb
    vira = orb.vira
    virb = orb.virb

    nCore = 100000000
    icore = np.empty(nCore, dtype=np.float64) # 800MB

    #------ get ab integral
    t1 = process_time()
    wt1 = time()
    # AA
    file = file_head + 'ABAA___'
    _write_1e_ao2mo_int(orb, ao_int, c_v_a, c_v_a, vira, vira, sym_A, icore, file, 1)
    if ref == 'UHF':
        # BB
        file = file_head + 'ABBB___'
        _write_1e_ao2mo_int(orb, ao_int, c_v_b, c_v_b, virb, virb, sym_A, icore, file, 1)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get VV integral:', t2-t1)
    print('wall time for get VV integral:', wt2-wt1)

    #------ get ij integral
    t1 = process_time()
    wt1 = time()
    # AA
    file = file_head + 'IJAA___'
    _write_1e_ao2mo_int(orb, ao_int, c_o_a, c_o_a, popa, popa, sym_A, icore, file, 1)
    if ref == 'UHF':
        # BB
        file = file_head + 'IJBB___'
        _write_1e_ao2mo_int(orb, ao_int, c_o_b, c_o_b, popb, popb, sym_A, icore, file, 1)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get OO integral:', t2-t1)
    print('wall time for get OO integral:', wt2-wt1)

    #------ get ia integral
    t1 = process_time()
    wt1 = time()
    # AA
    file = file_head + 'IAAA___'
    _write_1e_ao2mo_int(orb, ao_int, c_o_a, c_v_a, popa, vira, sym_A, icore, file, 2)
    if ref == 'UHF':
        # BB
        file = file_head + 'IABB___'
        _write_1e_ao2mo_int(orb, ao_int, c_o_b, c_v_b, popb, virb, sym_A, icore, file, 2)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get OV integral:', t2-t1)
    print('wall time for get OV integral:', wt2-wt1)

    #------ get ai integral
    t1 = process_time()
    wt1 = time()
    # AA
    file = file_head + 'AIAA___'
    _write_1e_ao2mo_int(orb, ao_int, c_v_a, c_o_a, vira, popa, sym_A, icore, file, 1)
    if ref == 'UHF':
        # BB
        file = file_head + 'AIBB___'
        _write_1e_ao2mo_int(orb, ao_int, c_v_b, c_o_b, virb, popb, sym_A, icore, file, 1)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get VO integral:', t2-t1)
    print('wall time for get VO integral:', wt2-wt1)

def _write_1e_ao2mo_int(orb, ao_int, c_p, c_q, pop_p, pop_q, sym_I, icore, file, type):
    '''
    type: 1 - ab, ij, ai; 2 - ia
    
    Calculate MO integral for specified r and s block
    '''

    # create file
    zorb = np.array([1,0,0,0,0,0,0,0])
    if type == 2:
        wi.creat_file(sym_I+1, pop_q, pop_p, 1, zorb, zorb, 1, file, icore)
    else:
        wi.creat_file(sym_I+1, pop_p, pop_q, 1, zorb, zorb, 1, file, icore)

    # mo(p,q) = ao(i,j)*c(ip)*c(jq)
    mo_int = lib.einsum('ij,ip,jq->pq', ao_int, c_p, c_q)
    print('norm of MO integral:',np.linalg.norm(mo_int))

    if type == 2:
        # ia->ai
        mo_int = mo_int.transpose(1,0)

    # change to Fortran 1D array
    mo_int = mo_int.ravel(order='F')

    # remove the zero block
    n_pq = wi.calc_num_of_pq(pop_p, pop_q, sym_I+1, 1)
    n_p = sum(pop_p)
    n_q = sum(pop_q)
    if type == 2:
        mo_int = wi.compress_by_symmerty(n_q, n_p, n_pq, pop_q, pop_p, sym_I+1, mo_int)
    else:
        mo_int = wi.compress_by_symmerty(n_p, n_q, n_pq, pop_p, pop_q, sym_I+1, mo_int)

    # store in file
    if type == 2:
        wi.write_block(n_pq, mo_int, sym_I+1, sym_I+1, pop_q, pop_p, 1, file, 1, icore)
    else:
        wi.write_block(n_pq, mo_int, sym_I+1, sym_I+1, pop_p, pop_q, 1, file, 1, icore)

def get_fock(orb, file_path=None):
    '''
    Write Fock matrix in file

    orb: orbital object
    '''
    if file_path is None:
        file_path = orb.tmpdir
    
    ref = orb.ref
    popa = orb.popa
    popb = orb.popb
    vira = orb.vira
    virb = orb.virb
    n_occ = orb.n_occ
    n_vir = orb.n_vir

    nCore = 100000000
    icore = np.empty(nCore, dtype=np.float64) # 800MB

    wi.init(orb.dirprd+1, orb.nirrep, file_path)

    #------ get ab integral
    t1 = process_time()
    wt1 = time()
    if ref == 'RHF':
        # AA
        file = 'FABAA___'
        fock = orb.fock[n_occ:,n_occ:]
        _write_1e_int(fock, vira, vira, icore, file, 1)
    elif ref == 'UHF':
        # AA
        file = 'FABAA___'
        fock = orb.fock[0][n_occ[0]:,n_occ[0]:]
        _write_1e_int(fock, vira, vira, icore, file, 1)
        # BB
        file = 'FABBB___'
        fock = orb.fock[1][n_occ[1]:,n_occ[1]:]
        _write_1e_int(fock, virb, virb, icore, file, 1)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get Fvv:', t2-t1)
    print('wall time for get Fvv:', wt2-wt1)

    #------ get ij integral
    t1 = process_time()
    wt1 = time()
    if ref == 'RHF':
        # AA
        file = 'FIJAA___'
        fock = orb.fock[:n_occ,:n_occ]
        _write_1e_int(fock, popa, popa, icore, file, 1)
    elif ref == 'UHF':
        # AA
        file = 'FIJAA___'
        fock = orb.fock[0][:n_occ[0],:n_occ[0]]
        _write_1e_int(fock, popa, popa, icore, file, 1)
        # BB
        file = 'FIJBB___'
        fock = orb.fock[1][:n_occ[1],:n_occ[1]]
        _write_1e_int(fock, popb, popb, icore, file, 1)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get Foo:', t2-t1)
    print('wall time for get Foo:', wt2-wt1)

    #------ get ia integral
    t1 = process_time()
    wt1 = time()
    if ref == 'RHF':
        # AA
        file = 'FAIAA___'
        fock = orb.fock[:n_occ,n_occ:]
        _write_1e_int(fock, popa, vira, icore, file, 2)
    elif ref == 'UHF':
        # AA
        file = 'FAIAA___'
        fock = orb.fock[0][:n_occ[0],n_occ[0]:]
        _write_1e_int(fock, popa, vira, icore, file, 2)
        # BB
        file = 'FAIBB___'
        fock = orb.fock[1][:n_occ[1],n_occ[1]:]
        _write_1e_int(fock, popb, virb, icore, file, 2)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get Fov:', t2-t1)
    print('wall time for get Fov:', wt2-wt1)

def _write_1e_int(mo_int, pop_p, pop_q, icore, file, type):
    '''
    Save 1e integral in DPD formate

    mo_int: MO integral without symmetry
    file: DPD file
    type: 1 - ab, ij; 2 - ia
    '''

    # create file
    zorb = np.array([1,0,0,0,0,0,0,0])
    if type == 2:
        wi.creat_file(1, pop_q, pop_p, 1, zorb, zorb, 1, file, icore)
    else:
        wi.creat_file(1, pop_p, pop_q, 1, zorb, zorb, 1, file, icore)

    if type == 2:
        # ia->ai
        mo_int = mo_int.transpose(1,0)

    # change to Fortran 1D array
    mo_int = mo_int.ravel(order='F')

    # remove the zero block
    n_pq = wi.calc_num_of_pq(pop_p, pop_q, 1, 1)
    n_p = sum(pop_p)
    n_q = sum(pop_q)
    if type == 2:
        mo_int = wi.compress_by_symmerty(n_q, n_p, n_pq, pop_q, pop_p, 1, mo_int)
    else:
        mo_int = wi.compress_by_symmerty(n_p, n_q, n_pq, pop_p, pop_q, 1, mo_int)

    # store in file
    if type == 2:
        wi.write_block(n_pq, mo_int, 1, 1, pop_q, pop_p, 1, file, 1, icore)
    else:
        wi.write_block(n_pq, mo_int, 1, 1, pop_p, pop_q, 1, file, 1, icore)

def _read_tmpdir():
    tmpdir='./tmp/'
    if Path('CCSOINP').exists():
        with open('CCSOINP', 'r') as f:
            lines = f.readlines()
        for line in lines:
            if line[0]=='!' or line[0]=='#':
               continue
            if 'tmpdir' in line:
                tmp_value = str(line)
                test_value=tmp_value.split('=')
                tmpdir = str(test_value[1]).strip()
    return tmpdir

def write_interface(orb):
    '''
    orb: orbital class
    '''
    # copy form pyscf_ccso/ccso_int and cc_pyscf/cc_int
    dirprd = orb.dirprd + 1
    inter.add_node_i2('dirprd', dirprd, 8, 8)

    if orb.ref == 'UHF':
        vrt = np.array(orb.vir).transpose(1,0)
        pop = np.array(orb.pop).transpose(1,0)
    else:
        vrt = np.array([orb.vir, orb.vir]).transpose(1,0)
        pop = np.array([orb.pop, orb.pop]).transpose(1,0)
    inter.add_node_i2('vrt', vrt, 8, 2)
    inter.add_node_i2('pop', pop, 8, 2)

    #isymmo = isymmo + 1
    #inter.add_node_i3('isymmo',isymmo,nmo,8,2)

    inter.add_node_iv('symrx', orb.symR[0]+1)
    inter.add_node_iv('symry', orb.symR[1]+1)
    inter.add_node_iv('symrz', orb.symR[2]+1)

    inter.add_node_iv('symx', orb.symX[0]+1)
    inter.add_node_iv('symy', orb.symX[1]+1)
    inter.add_node_iv('symz', orb.symX[2]+1)

    inter.add_node_dv('Escf', orb.e_scf)
    inter.add_node_str('group_n', orb.compGroup, len(orb.compGroup))
    inter.add_node_iv('frozen_orbit', orb.fz)
    inter.add_node_iv('nirrep', orb.nirrep)

    i = 0
    for key in orb.irrep_table.keys():
        inter.add_node_str(f'irrep_name{i}', key, len(key))
        i+=1

    if orb.ref == 'UHF':
        inter.add_node_d('orb_energy_alpha_sym', orb.mo_energy[0], len(orb.mo_energy[0]))
        inter.add_node_d('orb_energy_beta_sym', orb.mo_energy[1], len(orb.mo_energy[1]))
    else:
        inter.add_node_d('orb_energy_alpha_sym', orb.mo_energy, len(orb.mo_energy))
        inter.add_node_d('orb_energy_beta_sym', orb.mo_energy, len(orb.mo_energy))

    inter.add_node_iv('iprcln', 32768)
    inter.add_node_iv('iprcwd', 4096)
    inter.add_node_iv('ichcsz', 4096)

    print('finished save_ccso_int')

def get_mo_ovlp(orb, file_path=None):
    '''
    orb: orbital object
    
    calculate MO overlop matrix between alpha and beta orbital
    SoO_____: <i|J>
    SvV_____: <a|B>
    SvO_____: <a|I>
    SoV_____: <i|A>
    '''
    if file_path is None:
        file_path = orb.tmpdir

    t1 = process_time()
    wt1 = time()
    ref = orb.ref
    popa = orb.popa
    popb = orb.popb
    vira = orb.vira
    virb = orb.virb
    if ref == 'RHF':
        noa = orb.n_occ
        nob = orb.n_occ
        nva = orb.n_vir
        nvb = orb.n_vir
    elif ref == 'UHF':
        noa = orb.n_occ[0]
        nob = orb.n_occ[1]
        nva = orb.n_vir[0]
        nvb = orb.n_vir[1]

    # get AO integral <p|A|q>
    ao_ovlp = orb.mol.intor('int1e_ovlp')
    print('1e ao int shape',ao_ovlp.shape)
    if ao_ovlp.size != orb.norb_no_fz**2:
        print('ERROR: size of ao integral error', ao_ovlp.size, orb.norb_no_fz**2)

    # S_pQ = S_uv C_up C_vQ
    if ref == 'RHF':
        mo_ovlp = lib.einsum('uv,up,vq->pq',ao_ovlp,orb.mo_coeff,orb.mo_coeff)
    elif ref == 'UHF':
        mo_ovlp = lib.einsum('uv,up,vq->pq',ao_ovlp,orb.mo_coeff[0],orb.mo_coeff[1])

    # initialize mod
    wi.init(orb.dirprd+1, orb.nirrep, file_path)

    nCore = 100000000
    icore = np.empty(nCore, dtype=np.float64) # 800MB

    #------ S_iJ
    file = 'SoO_____'
    S = mo_ovlp[:noa,:nob]
    _write_1e_int(S, popa, popb, icore, file, 1)

    #------ S_aB
    file = 'SvV_____'
    S = mo_ovlp[noa:,nob:]
    _write_1e_int(S, vira, virb, icore, file, 1)

    #------ S_iA
    file = 'SoV_____'
    S = mo_ovlp[:noa,nob:]
    _write_1e_int(S, popa, virb, icore, file, 1)

    #------ S_aI
    file = 'SvO_____'
    S = mo_ovlp[noa:,:nob]
    _write_1e_int(S, vira, popb, icore, file, 1)

    t2 = process_time()
    wt2 = time()
    print('cpu time for get MO overlop matrix:', t2-t1)
    print('wall time for get MO overlop matrix:', wt2-wt1)

def get_mo_hcore(orb, file_path=None):
    '''
    orb: orbital object
    
    Write Fock matrix in file
    '''
    if file_path is None:
        file_path = orb.tmpdir

    ref = orb.ref
    popa = orb.popa
    popb = orb.popb
    vira = orb.vira
    virb = orb.virb
    n_occ = orb.n_occ
    n_vir = orb.n_vir

    nCore = 100000000
    icore = np.empty(nCore, dtype=np.float64) # 800MB

    wi.init(orb.dirprd+1, orb.nirrep, file_path)

    ao_hcore = scf.hf.get_hcore(orb.mol)
    if ref == 'RHF':
        mo_hcore = lib.einsum('uv,up,vq',ao_hcore,orb.mo_coeff,orb.mo_coeff)
    elif ref == 'UHF':
        mo_hcore = [None]*2
        mo_hcore[0] = lib.einsum('uv,up,vq',ao_hcore,orb.mo_coeff[0],orb.mo_coeff[0])
        mo_hcore[1] = lib.einsum('uv,up,vq',ao_hcore,orb.mo_coeff[1],orb.mo_coeff[1])

    #------ get ab integral
    t1 = process_time()
    wt1 = time()
    if ref == 'RHF':
        # AA
        file = 'hABAA___'
        fock = mo_hcore[n_occ:,n_occ:]
        _write_1e_int(fock, vira, vira, icore, file, 1)
    elif ref == 'UHF':
        # AA
        file = 'hABAA___'
        fock = mo_hcore[0][n_occ[0]:,n_occ[0]:]
        _write_1e_int(fock, vira, vira, icore, file, 1)
        # BB
        file = 'hABBB___'
        fock = mo_hcore[1][n_occ[1]:,n_occ[1]:]
        _write_1e_int(fock, virb, virb, icore, file, 1)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get hvv:', t2-t1)
    print('wall time for get hvv:', wt2-wt1)

    #------ get ij integral
    t1 = process_time()
    wt1 = time()
    if ref == 'RHF':
        # AA
        file = 'hIJAA___'
        fock = mo_hcore[:n_occ,:n_occ]
        _write_1e_int(fock, popa, popa, icore, file, 1)
    elif ref == 'UHF':
        # AA
        file = 'hIJAA___'
        fock = mo_hcore[0][:n_occ[0],:n_occ[0]]
        _write_1e_int(fock, popa, popa, icore, file, 1)
        # BB
        file = 'hIJBB___'
        fock = mo_hcore[1][:n_occ[1],:n_occ[1]]
        _write_1e_int(fock, popb, popb, icore, file, 1)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get hoo:', t2-t1)
    print('wall time for get hoo:', wt2-wt1)

    #------ get ia integral
    t1 = process_time()
    wt1 = time()
    if ref == 'RHF':
        # AA
        file = 'hAIAA___'
        fock = mo_hcore[:n_occ,n_occ:]
        _write_1e_int(fock, popa, vira, icore, file, 2)
    elif ref == 'UHF':
        # AA
        file = 'hAIAA___'
        fock = mo_hcore[0][:n_occ[0],n_occ[0]:]
        _write_1e_int(fock, popa, vira, icore, file, 2)
        # BB
        file = 'hAIBB___'
        fock = mo_hcore[1][:n_occ[1],n_occ[1]:]
        _write_1e_int(fock, popb, virb, icore, file, 2)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get hov:', t2-t1)
    print('wall time for get hov:', wt2-wt1)

def get_2e_int_outcore(orb, file_path=None, integral='int2e', file_head=' ', sym_A=0):
    '''
    orb: orbital object
    integral: str, the name of integral need get, see pyscf/gto/moleintor.py
    sym_A: symmerty of operator A
    file_head: str with len=1

    Get 2e integral A^pq_rs = <pq|A|rs>-<pq|A|sr> and write in hard disk in DPD formate.
    PySCF store 2e-integral as (pq|rs)=<pr|qs>
    '''
    if file_path is None:
        file_path = orb.tmpdir

    mol = orb.mol
    ref = orb.ref

    # get orbital coefficient
    if ref == 'RHF':
        c_o_a = orb.mo_coeff[:, 0:orb.n_occ]
        c_v_a = orb.mo_coeff[:, orb.n_occ:]
        c_o_b = c_o_a
        c_v_b = c_v_a
        n_o_a = orb.n_occ
        n_o_b = orb.n_occ
        n_v_a = orb.n_vir
        n_v_b = orb.n_vir
    elif ref == 'UHF':
        c_o_a = orb.mo_coeff[0][:, 0:orb.n_occ[0]]
        c_v_a = orb.mo_coeff[0][:, orb.n_occ[0]:]
        c_o_b = orb.mo_coeff[1][:, 0:orb.n_occ[1]]
        c_v_b = orb.mo_coeff[1][:, orb.n_occ[1]:]
        n_o_a = orb.n_occ[0]
        n_o_b = orb.n_occ[1]
        n_v_a = orb.n_vir[0]
        n_v_b = orb.n_vir[1]
    assert c_o_a.shape[1] == n_o_a
    print('shape of MO coeff occupied alpha:', c_o_a.shape)
    assert c_o_b.shape[1] == n_o_b
    print('shape of MO coeff virtual alpha:', c_v_a.shape)
    assert c_v_a.shape[1] == n_v_a
    print('shape of MO coeff occupied beta:', c_o_b.shape)
    assert c_v_b.shape[1] == n_v_b
    print('shape of MO coeff virtual beta:', c_v_b.shape)

    #------ get 2e integral file
    h5file = file_path+'2e.dat'
    print('size of va^4 (GiB):', n_v_a**4*8/1024**3)
    print('size of vb^4 (GiB):', n_v_b**4*8/1024**3)
    print('size of va^3*oa (GiB):', n_v_a**3*n_o_a*8/1024**3)
    print('size of vb^3*ob (GiB):', n_v_b**3*n_o_b*8/1024**3)
    print('size of va^2*oa^2 (GiB):', n_v_a**2*n_o_a**2*8/1024**3)
    print('size of vb^2*ob^2 (GiB):', n_v_b**2*n_o_b**2*8/1024**3)

    t1 = process_time()
    wt1 = time()
    #----------------------------------------------------
    # <pq|rs> need cycle for s, and in I/O procedure,
    # it is best to cycle for the first idx u in (uv,wx)
    # so we store <pq|rs> as (sq,rp):
    #   <pq|rs> = (pr|qs) = (qs|pr) = (sq|rp)
    # the final equation need the s4 symmerty,
    # this symmerty require operator A be Hermitian operator
    #----------------------------------------------------
    print('ATTENTATION: the 2e integral need satisfy s4 symmetry - (pq|rs)=(qp|rs)=(pq|sr)=(qp|sr)')
    #--- <ab||cd>
    # <aB|cD> = (ac|BD) = (BD|ac) = (DB|ca)
    ao2mo.outcore.general(mol, (c_v_b,c_v_b,c_v_a,c_v_a),intor=integral,\
                          erifile=h5file, dataname='vVvV', compact=False)
    if ref == 'UHF':
        # <ab|cd> = (ac|bd) = (bd|ac) = (db|ca)
        ao2mo.outcore.general(mol, (c_v_a,c_v_a,c_v_a,c_v_a),intor=integral,\
                              erifile=h5file, dataname='vvvv', compact=False)
        # <AB|CD> = (AC|BD) = (BD|AC) = (DB|CA)
        ao2mo.outcore.general(mol, (c_v_b,c_v_b,c_v_b,c_v_b),intor=integral,\
                              erifile=h5file, dataname='VVVV', compact=False)
    #--- <ij||kl>
    # <iJ|kL> = (ik|JL) = (JL|ik) = (LJ|ki)
    ao2mo.outcore.general(mol, (c_o_b,c_o_b,c_o_a,c_o_a),intor=integral,\
                          erifile=h5file, dataname='oOoO', compact=False)
    if ref == 'UHF':
        # <ij|kl> = (ik|jl) = (jl|ik) = (lj|ki)
        ao2mo.outcore.general(mol, (c_o_a,c_o_a,c_o_a,c_o_a),intor=integral,\
                              erifile=h5file, dataname='oooo', compact=False)
        # <IJ|KL> = (IK|JL) = (JL|IK) = (LJ|KI)
        ao2mo.outcore.general(mol, (c_o_b,c_o_b,c_o_b,c_o_b),intor=integral,\
                              erifile=h5file, dataname='OOOO', compact=False)
    #--- <ab||ij>
    # <aB|iJ> = (ai|BJ) = (BJ|ai) = (JB|ia)
    ao2mo.outcore.general(mol, (c_o_b,c_v_b,c_o_a,c_v_a),intor=integral,\
                          erifile=h5file, dataname='vVoO', compact=False)
    if ref == 'UHF':
        # <ab|ij> = (ai|bj) = (bj|ai) = (jb|ia)
        ao2mo.outcore.general(mol, (c_o_a,c_v_a,c_o_a,c_v_a),intor=integral,\
                              erifile=h5file, dataname='vvoo', compact=False)
        # <AB|IJ> = (AI|BJ) = (BJ|AI) = (JB|IA)
        ao2mo.outcore.general(mol, (c_o_b,c_v_b,c_o_b,c_v_b),intor=integral,\
                              erifile=h5file, dataname='VVOO', compact=False)
    #--- <ab||ci>
    # <aB|cI> = (ac|BI) = (BI|ac) = (IB|ca)
    ao2mo.outcore.general(mol, (c_o_b,c_v_b,c_v_a,c_v_a),intor=integral,\
                          erifile=h5file, dataname='vVvO', compact=False)
    if ref == 'UHF':
        # <ab|ci> = (ac|bi) = (bi|ac) = (ib|ca)
        ao2mo.outcore.general(mol, (c_o_a,c_v_a,c_v_a,c_v_a),intor=integral,\
                              erifile=h5file, dataname='vvvo', compact=False)
        # <AB|CI> = (AC|BI) = (BI|AC) = (IB|CA)
        ao2mo.outcore.general(mol, (c_o_b,c_v_b,c_v_b,c_v_b),intor=integral,\
                              erifile=h5file, dataname='VVVO', compact=False)
        # <Ab|Ci> = (AC|bi) = (bi|AC) = (ib|CA)
        ao2mo.outcore.general(mol, (c_o_a,c_v_a,c_v_b,c_v_b),intor=integral,\
                              erifile=h5file, dataname='VvVo', compact=False)
    #--- <ij||ka>
    # <iJ|kA> = (ik|JA) = (JA|ik) = (AJ|ki)
    ao2mo.outcore.general(mol, (c_v_b,c_o_b,c_o_a,c_o_a),intor=integral,\
                          erifile=h5file, dataname='oOoV', compact=False)
    if ref =='UHF':
        # <ij|ka> = (ik|ja) = (ja|ik) = (aj|ki)
        ao2mo.outcore.general(mol, (c_v_a,c_o_a,c_o_a,c_o_a),intor=integral,\
                              erifile=h5file, dataname='ooov', compact=False)
        # <IJ|KA> = (IK|JA) = (JA|IK) = (AJ|KI)
        ao2mo.outcore.general(mol, (c_v_b,c_o_b,c_o_b,c_o_b),intor=integral,\
                              erifile=h5file, dataname='OOOV', compact=False)
        # <Ij|Ka> = (IK|ja) = (ja|IK) = (aj|KI)
        ao2mo.outcore.general(mol, (c_v_a,c_o_a,c_o_b,c_o_b),intor=integral,\
                              erifile=h5file, dataname='OoOv', compact=False)
    #--- <ia||bj>
    # <ia|bj> = (ib|aj) = (aj|ib) = (ja|bi)
    ao2mo.outcore.general(mol, (c_o_a,c_v_a,c_v_a,c_o_a),intor=integral,\
                          erifile=h5file, dataname='ovvo', compact=False)
    # <ai|bj> = (ab|ij) = (ij|ab) = (ji|ba)
    ao2mo.outcore.general(mol, (c_o_a,c_o_a,c_v_a,c_v_a),intor=integral,\
                          erifile=h5file, dataname='vovo', compact=False)
    # <iA|bJ> = (ib|AJ) = (AJ|ib) = (JA|bi)
    ao2mo.outcore.general(mol, (c_o_b,c_v_b,c_v_a,c_o_a),intor=integral,\
                          erifile=h5file, dataname='oVvO', compact=False)
    # <aI|bJ> = (ab|IJ) = (IJ|ab) = (JI|ba)
    ao2mo.outcore.general(mol, (c_o_b,c_o_b,c_v_a,c_v_a),intor=integral,\
                          erifile=h5file, dataname='vOvO', compact=False)
    if ref == 'UHF':
        # <IA|BJ> = (IB|AJ) = (AJ|IB) = (JA|BI)
        ao2mo.outcore.general(mol, (c_o_b,c_v_b,c_v_b,c_o_b),intor=integral,\
                              erifile=h5file, dataname='OVVO', compact=False)
        # <AI|BJ> = (AB|IJ) = (IJ|AB) = (JI|BA)
        ao2mo.outcore.general(mol, (c_o_b,c_o_b,c_v_b,c_v_b),intor=integral,\
                              erifile=h5file, dataname='VOVO', compact=False)
        # <Ai|Bj> = (AB|ij) = (ij|AB) = (ji|BA)
        ao2mo.outcore.general(mol, (c_o_a,c_o_a,c_v_b,c_v_b),intor=integral,\
                              erifile=h5file, dataname='VoVo', compact=False)
    t2 = process_time()
    wt2 = time()
    print('cpu time for 2e int AO->MO transform:', t2-t1)
    print('wall time for 2e int AO->MO transform:', wt2-wt1)

    # initialize mod
    wi.init(orb.dirprd+1, orb.nirrep, file_path)

    f = h5py.File(h5file,'r')
    _2e_int_outcore(orb, f, file_head, sym_A)

def _2e_int_outcore(orb, h5obj, file_head=' ', sym_A=0):
    '''
    h5obj: hdf5 object, h5obj = h5py.File(h5file,'r')
    sym_A: symmetry of 2e integral
    '''
    ref = orb.ref
    popa = orb.popa
    popb = orb.popb
    vira = orb.vira
    virb = orb.virb

    nCore = 100000000
    icore = np.empty(nCore, dtype=np.float64) # 800MB

    if file_head == ' ':
        default_name = True
    else:
        default_name = False

    #------ get abcd integral
    t1 = process_time()
    wt1 = time()
    # <aB||cD> = <aB|cD>
    if default_name:
        file = 'ABCDAB__'
    else:
        file = file_head + '2e1ABAB'
    _pre_pqr_block(orb, vira, virb, 1, vira, virb, 1, \
                   h5obj, 'vVvV', 2, False, sym_A, icore, file)
    if ref == 'UHF':
        # <ab||cd> = <ab|cd>-<ba|cd>
        if default_name:
            file = 'ABCDAA__'
        else:
            file = file_head + '2e1AAAA'
        _pre_pqr_block(orb, vira, vira, 0, vira, vira, 0, \
                       h5obj, 'vvvv', 1, False, sym_A, icore, file)
        # <AB||CD> = <AB|CD>-<BA|CD>
        if default_name:
            file = 'ABCDBB__'
        else:
            file = file_head + '2e1BBBB'
        _pre_pqr_block(orb, virb, virb, 0, virb, virb, 0, \
                       h5obj, 'VVVV', 1, False, sym_A, icore, file)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get VVVV integral:', t2-t1)
    print('wall time for get VVVV integral:', wt2-wt1)

    #------ get abij integral
    t1 = process_time()
    wt1 = time()
    # <aB||iJ> = <aB|iJ>
    if default_name:
        file = 'ABIJAB__'
    else:
        file = file_head + '2e2ABAB'
    _pre_pqr_block(orb, vira, virb, 1, popa, popb, 1, \
                   h5obj, 'vVoO', 2, False, sym_A, icore, file)
    if ref == 'UHF':
        # <ab||ij> = <ab|ij>-<ba|ij>
        if default_name:
            file = 'ABIJAA__'
        else:
            file = file_head + '2e2AAAA'
        _pre_pqr_block(orb, vira, vira, 0, popa, popa, 0, \
                       h5obj, 'vvoo', 1, False, sym_A, icore, file)
        # <AB||IJ> = <AB|IJ>-<BA|IJ>
        if default_name:
            file = 'ABIJBB__'
        else:
            file = file_head + '2e2BBBB'
        _pre_pqr_block(orb, virb, virb, 0, popb, popb, 0, \
                       h5obj, 'VVOO', 1, False, sym_A, icore, file)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get VVOO integral:', t2-t1)
    print('wall time for get VVOO integral:', wt2-wt1)

    #------ get ijkl integral
    t1 = process_time()
    wt1 = time()
    # <iJ||kL>=<iJ|kL>
    if default_name:
        file = 'IJKLAB__'
    else:
        file = file_head + '2e3ABAB'
    _pre_pqr_block(orb, popa, popb, 1, popa, popb, 1, \
                   h5obj, 'oOoO', 2, False, sym_A, icore, file)
    if ref == 'UHF':
        # <ij||kl>=<ij|kl>-<ji|kl>
        if default_name:
            file = 'IJKLAA__'
        else:
            file = file_head + '2e3AAAA'
        _pre_pqr_block(orb, popa, popa, 0, popa, popa, 0, \
                       h5obj, 'oooo', 1, False, sym_A, icore, file)
        # <IJ||KL>=<IJ|KL>-<JI|KL>
        if default_name:
            file = 'IJKLBB__'
        else:
            file = file_head + '2e3BBBB'
        _pre_pqr_block(orb, popb, popb, 0, popb, popb, 0, \
                       h5obj, 'OOOO', 1, False, sym_A, icore, file)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get OOOO integral:', t2-t1)
    print('wall time for get OOOO integral:', wt2-wt1)

    #------ get abci integral
    t1 = process_time()
    wt1 = time()
    # <aB||cI> = <aB|cI>
    if default_name:
        file = 'ABCIAB__'
    else:
        file = file_head + '2e4ABAB'
    _pre_pqr_block(orb, vira, virb, 1, vira, popb, 1, \
                   h5obj, 'vVvO', 2, False, sym_A, icore, file)
    if ref == 'UHF':
        # <ab||ci> = <ab|ci>-<ba|ci>
        if default_name:
            file = 'ABCIAA__'
        else:
            file = file_head + '2e4AAAA'
        _pre_pqr_block(orb, vira, vira, 0, vira, popa, 1, \
                       h5obj, 'vvvo', 1, False, sym_A, icore, file)
        # <AB||CI> = <AB|CI>-<BA|CI>
        if default_name:
            file = 'ABCIBB__'
        else:
            file = file_head + '2e4BBBB'
        _pre_pqr_block(orb, virb, virb, 0, virb, popb, 1, \
                       h5obj, 'VVVO', 1, False, sym_A, icore, file)
        # <Ab||Ci> = <Ab|Ci>
        if default_name:
            file = 'ABCIBA__'
        else:
            file = file_head + '2e4BABA'
        _pre_pqr_block(orb, virb, vira, 1, virb, popa, 1, \
                       h5obj, 'VvVo', 2, False, sym_A, icore, file)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get VVVO integral:', t2-t1)
    print('wall time for get VVVO integral:', wt2-wt1)

    #------ get ijka integral
    t1 = process_time()
    wt1 = time()
    # <iJ||kA> = <iJ|kA>
    if default_name:
        file = 'IJKAAB__'
    else:
        file = file_head + '2e5ABAB'
    _pre_pqr_block(orb, popa, popb, 1, popa, virb, 1, \
                   h5obj, 'oOoV', 2, False, sym_A, icore, file)
    if ref == 'UHF':
        # <ij||ka> = <ij|ka>-<ji|ka>
        if default_name:
            file = 'IJKAAA__'
        else:
            file = file_head + '2e5AAAA'
        _pre_pqr_block(orb, popa, popa, 0, popa, vira, 1, \
                       h5obj, 'ooov', 1, False, sym_A, icore, file)
        # <IJ||KA> = <IJ|KA>-<JI|KA>
        if default_name:
            file = 'IJKABB__'
        else:
            file = file_head + '2e5BBBB'
        _pre_pqr_block(orb, popb, popb, 0, popb, virb, 1, \
                       h5obj, 'OOOV', 1, False, sym_A, icore, file)
        # <Ij||Ka> = <Ij|Ka>
        if default_name:
            file = 'IJKABA__'
        else:
            file = file_head + '2e5BABA'
        _pre_pqr_block(orb, popb, popa, 1, popb, vira, 1, \
                       h5obj, 'OoOv', 2, False, sym_A, icore, file)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get OOOV integral:', t2-t1)
    print('wall time for get OOOV integral:', wt2-wt1)

    #------ get iabj integral
    t1 = process_time()
    wt1 = time()
    # <ia|bj> = <ia|bj>-<ai|bj>
    if default_name:
        file = 'BIAJAAAA'
    else:
        file = file_head + '2e6AAAA'
    _pre_pqr_block(orb, popa, vira, 1, vira, popa, 1, \
                   h5obj, ['ovvo','vovo'], 3, True, sym_A, icore, file)
    # <iA||bJ> = <iA|bJ>
    if default_name:
        file = 'BIAJAABB'
    else:
        file = file_head + '2e6AABB'
    _pre_pqr_block(orb, popa, virb, 1, vira, popb, 1, \
                   h5obj, 'oVvO', 2, True, sym_A, icore, file)
    # <Ia||bJ> = -<aI|bJ>
    if default_name:
        file = 'BIAJABAB'
    else:
        file = file_head + '2e6ABAB'
    _pre_pqr_block(orb, popb, vira, 1, vira, popb, 1, \
                   h5obj, 'vOvO', 4, True, sym_A, icore, file)
    if ref == 'UHF':
        # <IA|BJ> = <IA|BJ>-<AI|BJ>
        if default_name:
            file = 'BIAJBBBB'
        else:
            file = file_head + '2e6BBBB'
        _pre_pqr_block(orb, popb, virb, 1, virb, popb, 1, \
                       h5obj, ['OVVO','VOVO'], 3, True, sym_A, icore, file)
        # <iA||Bj> = -<Ai|Bj>
        if default_name:
            file = 'BIAJBABA'
        else:
            file = file_head + '2e6BABA'
        _pre_pqr_block(orb, popa, virb, 1, virb, popa, 1, \
                       h5obj, 'VoVo', 4, True, sym_A, icore, file)
    t2 = process_time()
    wt2 = time()
    print('cpu time for get OVVO integral:', t2-t1)
    print('wall time for get OVVO integral:', wt2-wt1)

def _pre_pqr_block(orb, pop_p, pop_q, ioppq, pop_r, pop_s, ioprs,\
                   h5obj, int_type, type, iabj, sym_I, icore, file):
    '''
    orb: orbital object
    h5obj: hdf5 object, h5obj = h5py.File(h5file,'r')
    int_type: str, dataname for h5file; for type=3 need be a [int_pq, int_qp]
    type: int, 1 - <pq||rs> = <pq|rs>-<qp|rs>, pq with same type
               2 - <pq||rs> = <pq|rs>
               3 - <pq||rs> = <pq|rs>-<qp|rs>, pq with different type
               4 - <pq||rs> = -<qp|rs>
    iabj: bool, whether be iabj integral, need resort as (bi,aj)
    sym_I: int, symmerty of <pq||rs>
    file: str, DPD format file
    '''
    nirrep = orb.nirrep
    n_p = sum(pop_p[0:nirrep])
    n_q = sum(pop_q[0:nirrep])
    n_r = sum(pop_r[0:nirrep])
    n_s = sum(pop_s[0:nirrep])
    if iabj:
        wi.creat_file(sym_I+1, pop_r, pop_p, 1, pop_q, pop_s, 1, file, icore)
    else:
        wi.creat_file(sym_I+1, pop_p, pop_q, ioppq, pop_r, pop_s, ioprs, file, icore)

    #------ cycle for sym_s
    for sym_s in range(0, nirrep):

        #------ cycle for symmetry adapted s
        for s in range(0, pop_s[sym_s]):
            s_ns = sum(pop_s[0:sym_s])+s # s without symmerty

            #--- get for a s
            if type in [1,2]: 
                # 1 - <pq||rs> = <pq|rs>-<qp|rs>, pq with same type
                # 2 - <pq||rs> = <pq|rs>
                # get <pq|rs> = (sq|rp)(*,q,r,p)
                block = _h5_read_qrs(h5obj,int_type,s_ns,n_s,n_q,n_r,n_p)
                # <pq|rs>(q,r,p) -> <pq|rs>(r,p,q)
                block = np.ascontiguousarray(block.transpose(1,2,0))
                if type == 1:
                    # <pq|rs>(r,p,q)=<qp|rs>(r,q,p)-><qp|rs>(r,p,q)
                    block2 = np.ascontiguousarray(block.transpose(0,2,1))
                    # <pq||rs>(r,p,q) = <pq|rs>(r,p,q)-<qp|rs>(r,p,q)
                    block = block-block2
            elif type == 3: # <pq||rs> = <pq|rs>-<qp|rs>, pq with different type
                # get <pq|rs>=(sq|rp)(*,q,r,p)
                block = _h5_read_qrs(h5obj,int_type[0],s_ns,n_s,n_q,n_r,n_p)
                # <pq|rs>(q,r,p) -> <pq|rs>(r,p,q)
                block = np.ascontiguousarray(block.transpose(1,2,0))
                # get <qp|rs>=(sp|rq)(*,p,r,q)
                block2 = _h5_read_qrs(h5obj,int_type[1],s_ns,n_s,n_p,n_r,n_q)
                # <qp|rs>(p,r,q) -> <qp|rs>(r,p,q)
                block2 = np.ascontiguousarray(block2.transpose(1,0,2))
                # <pq||rs>(r,p,q) = <pq|rs>(r,p,q)-<qp|rs>(r,p,q)
                block = block - block2
            elif type == 4: # <pq||rs> = -<qp|rs>
                # get <qp|rs>=(sp|rq)(p,r,q)
                block = _h5_read_qrs(h5obj,int_type,s_ns,n_s,n_p,n_r,n_q)
                # <qp|rs>(p,r,q) -> <qp|rs>(r,p,q)
                # <pq||rs>(r,p,q) = -<qp|rs>(r,p,q)
                block = -np.ascontiguousarray(block.transpose(1,0,2))

            if iabj:
                # iabj(b,i,a)->biaj(a,b,i)
                block = np.ascontiguousarray(block.transpose(2,0,1))
                _write_pqr_block(orb, block, pop_r, pop_p, 1, pop_q, pop_s, 1,\
                           sym_s, s_ns, sym_I, icore, file)
            else:
                _write_pqr_block(orb, block, pop_p, pop_q, ioppq, pop_r, pop_s, ioprs,\
                           sym_s, s_ns, sym_I, icore, file)

def _write_pqr_block(orb, mo_int, pop_p, pop_q, ioppq, pop_r, pop_s, ioprs,\
               sym_s, s, sym_I, icore, file):
    '''
    mo_int: 3D array, I(pq,r*),sort as (r,pq)
    s: s without symmetry
    file: DPD format file
    
    write 2e int I(pq,r*) for a specified s
    '''
    nirrep = orb.nirrep
    dirprd = orb.dirprd
    s_sa = s-sum(pop_s[0:sym_s]) # symmetry adapted s

    #------ cycle for sym_r
    for sym_r in range(0,nirrep):
        sym_rs = dirprd[sym_r,sym_s]
        sym_pq = dirprd[sym_rs,sym_I]
        num_r = pop_r[sym_r]

        if ioprs == 1:
            num_r_sa = num_r
        elif ioprs == 0:
            if sym_r < sym_s:
                num_r_sa = num_r
            elif sym_r == sym_s:
                num_r_sa = s_sa
            else:
                continue
        else:
            print('ERROR: ioprs error:',ioprs)
            exit()

        #------ cycle for symmerty adapted r
        for r_sa in range(0,num_r_sa):
            r = sum(pop_r[0:sym_r])+r_sa
            block = mo_int[r,:,:]
            rs = wi.put_rs(sym_rs+1, pop_r, pop_s, sym_r+1, sym_s+1, ioprs, r_sa+1, s_sa+1)
            rs -= 1
            _write_pq_block(orb, block, pop_p, pop_q, ioppq, sym_pq, sym_I, rs, icore, file)

def _write_pq_block(orb, mo_int, pop_p, pop_q, ioppq, sym_pq, sym_I, rs, icore, file):
    '''
    mo_int: 2D array, I(pq,**)
    pop_p/q: 1D array, symmertry population of p/q
    rs: symmetry adapted rs in I(pq,rs)
    sym_pq: symmetry of pq in I(pq,rs)
    sym_I: symmerty of I(pq,rs)
    file: DPD format file
    
    write 2e int I(pq,**)
    '''
    # change to Fortran 1D array
    mo_int = mo_int.ravel(order='F')
    
    # remove the zero block
    n_pq = wi.calc_num_of_pq(pop_p, pop_q, sym_pq+1, 1)
    n_p = sum(pop_p[0:orb.nirrep])
    n_q = sum(pop_q[0:orb.nirrep])
    mo_int = wi.compress_by_symmerty(n_p, n_q, n_pq, pop_p, pop_q, sym_pq+1, mo_int)
    
    # store in file
    wi.write_block(n_pq, mo_int, sym_pq+1, sym_I+1, pop_p, pop_q, ioppq, file, rs+1, icore)

def _h5_read_qrs(h5obj, datename, i1_ns, n1, n2, n3, n4):
    '''
    the chunks of (pq,rs) is (min(p,q),min(r,s))
    and q in pq, s in rs is continue
    so it is best to read for a sepcifed p or r
    but this need s4 symmetry: (pq|rs)=(qp|rs)=(pq|sr)=(qp|sr)

    h5obj: hdf5 object, h5obj = h5py.File(h5file,'r')
    i1_ns: p without symmerty in (pq|rs)
    n1: total number of p in (pq|rs)
    n2: total number of q in (pq|rs)
    n3: total number of r in (pq|rs)
    n4: total number of s in (pq|rs)
    
    return:
        (*q|rs)(q,r,s)
    '''
    start = i1_ns*n2
    end = start+n2
    qrs = h5obj[datename][start:end,:]

    return qrs.reshape(n2,n3,n4)

def scf_copy(mf):
    '''
    Create a new scf object which have the same type, mo_oceff, and mo_occ with mf

    mf: scf object
    '''
    mf_out = mf.copy()

    cls_name = str(type(mf))
    if 'UHF' in cls_name:
        mf_out.mo_coeff = [mf.mo_coeff[0].copy(), mf.mo_coeff[1].copy()]
        mf_out.mo_occ = [mf.mo_occ[0].copy(), mf.mo_occ[1].copy()]
    else:
        mf_out.mo_coeff = mf.mo_coeff.copy()
        mf_out.mo_occ = mf.mo_occ.copy()

    return mf_out

def get_pair_orbital_idx(orb_i, orb_f, ao_ovlp=None):
    '''
    Determine idx of orbital both cooupied by initial and final state
    and idx of excited core and target beta orbital.
    Save pocc_a, pocc_b, pvir_a, pvir_b, u_b, x_b in interface
    u is core orbital, x is target orbital

    orb_i: initial state, need be RHF state
    orb_f: final state, need be UHF state
    ao_ovlp: overlop matrix in AO basis
    '''

    if ao_ovlp is None:
        ao_ovlp = orb_i.mol.intor('int1e_ovlp')

    no_i = orb_i.n_occ
    nv_i = orb_i.n_vir
    noa_f = orb_f.n_occ[0]
    nob_f = orb_f.n_occ[1]
    nva_f = orb_f.n_vir[0]
    nvb_f = orb_f.n_vir[1]
    assert no_i==noa_f==nob_f
    assert nv_i==nva_f==nvb_f

    co_i = orb_i.mo_coeff[:,:no_i]
    cv_i = orb_i.mo_coeff[:,no_i:]
    coa_f = orb_f.mo_coeff[0][:,:noa_f]
    cva_f = orb_f.mo_coeff[0][:,noa_f:]
    cob_f = orb_f.mo_coeff[1][:,:nob_f]
    cvb_f = orb_f.mo_coeff[1][:,nob_f:]

    pocc_a = np.zeros((2,no_i),dtype=np.int64,order='F')
    pvir_a = np.zeros((2,nv_i),dtype=np.int64,order='F')
    pocc_b = np.zeros((2,no_i-1),dtype=np.int64,order='F')
    pvir_b = np.zeros((2,nv_i-1),dtype=np.int64,order='F')
    u_b = np.zeros(2,dtype=np.int64)
    x_b = np.zeros(2,dtype=np.int64)

    # <io|fo>alpha
    S = np.einsum('uv,up,vq->pq', ao_ovlp, co_i, coa_f)
    r,c = np.where(np.abs(S)>0.9999)
    assert r.shape[0]==no_i
    for i in range(0,no_i):
        pocc_a[0,i] = r[i]+1
        pocc_a[1,i] = c[i]+1
    print('pocc_a:\n',pocc_a)

    # <iv|fv>alpha
    S = np.einsum('uv,up,vq->pq', ao_ovlp, cv_i, cva_f)
    r,c = np.where(np.abs(S)>0.9999)
    assert r.shape[0]==nv_i
    for i in range(0,nv_i):
        pvir_a[0,i] = r[i]+1
        pvir_a[1,i] = c[i]+1
    print('pvir_a:\n',pvir_a)

    # <io|fo>beta
    S = np.einsum('uv,up,vq->pq', ao_ovlp, co_i, cob_f)
    r,c = np.where(np.abs(S)>0.9999)
    assert r.shape[0]==no_i-1
    for i in range(0,no_i-1):
        pocc_b[0,i] = r[i]+1
        pocc_b[1,i] = c[i]+1
    print('pocc_b:\n',pocc_b)

    # <iv|fv>beta
    S = np.einsum('uv,up,vq->pq', ao_ovlp, cv_i, cvb_f)
    r,c = np.where(np.abs(S)>0.9999)
    assert r.shape[0]==nv_i-1
    for i in range(0,nv_i-1):
        pvir_b[0,i] = r[i]+1
        pvir_b[1,i] = c[i]+1
    print('pvir_b:\n',pvir_b)

    # <io_beta|fv_beta>
    S = np.einsum('uv,up,vq->pq', ao_ovlp, co_i, cvb_f)
    r,c = np.where(np.abs(S)>0.9999)
    assert r.shape[0]==1
    u_b[0] = r[0]+1
    u_b[1] = c[0]+1
    print('u_b:',u_b)

    # <iv_beta|fo_beta>
    S = np.einsum('uv,up,vq->pq', ao_ovlp, cv_i, cob_f)
    r,c = np.where(np.abs(S)>0.9999)
    assert r.shape[0]==1
    x_b[0] = r[0]+1
    x_b[1] = c[0]+1
    print('x_b:',x_b)

    inter.add_node_i2('pocc_a', pocc_a, 2, no_i)
    inter.add_node_i2('pvir_a', pvir_a, 2, nv_i)
    inter.add_node_i2('pocc_b', pocc_b, 2, no_i-1)
    inter.add_node_i2('pvir_b', pvir_b, 2, nv_i-1)
    inter.add_node_i('u_b', u_b, 2)
    inter.add_node_i('x_b', x_b, 2)

def max_overlop_MO(current, target, ovlp):
    '''
    Rotate current orbital to make it have max overlop to the target orbital
    target: 2D-array, target orbital
    current: 2D-array, current orbital
    ovlp: AO overlop matirx

    return: mo_coeff after rotation
    '''

    # I minimize the ||(target)-(current)*C||_F use Procrustes algorithm
    # ref: G. H. Gloub, Matrix Computation, 4th, P329

    # S = <current|target>
    S_mo = np.einsum('uv,up,vq->pq', ovlp, current, target)
    print('Before rotation:\n',[S_mo[i,i] for i in range(0,S_mo.shape[0])])
    # SVD
    U,S,VT = np.linalg.svd(S_mo)
    # C = U*V^T
    C = np.einsum('ik,kj->ij', U, VT)
    # new_mo(u,i) = mo(u,p)*C(p,i)
    coeff = np.einsum('up,pi->ui', current, C)

    S_mo = np.einsum('uv,up,vq->pq', ovlp, coeff, target)
    print('After rotation:\n',[S_mo[i,i] for i in range(0,S_mo.shape[0])])

    return coeff

def get_mass_center(mol):
    '''
    Calculate mass center for a moleculer.
    '''
    mass = mol.atom_mass_list()
    coords = mol.atom_coords()
    mass_center = np.einsum('i,ij->j', mass, coords)/mass.sum()
    print(f'Mass center (a.u.): {mass_center[0]:.8f} {mass_center[1]:.8f} {mass_center[2]:.8f}')

    return mass_center
