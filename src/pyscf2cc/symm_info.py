import numpy as np
from pyscf import scf, symm

# direct prduct table used in PySCF
dirprd = np.array(((0, 1, 2, 3, 4, 5, 6, 7),
                   (1, 0, 3, 2, 5, 4, 7, 6),
                   (2, 3, 0, 1, 6, 7, 4, 5),
                   (3, 2, 1, 0, 7, 6, 5, 4),
                   (4, 5, 6, 7, 0, 1, 2, 3),
                   (5, 4, 7, 6, 1, 0, 3, 2),
                   (6, 7, 4, 5, 2, 3, 0, 1),
                   (7, 6, 5, 4, 3, 2, 1, 0)),
                   dtype = np.int64)

def canonical(mf):
    '''
    mf: scf object

    get canonical orbital
    this function will change mf object
    '''
    energy, coeff = mf.canonicalize(mf.mo_coeff, mf.mo_occ)
    mf.mo_energy = energy
    mf.mo_coeff = coeff
    print('orb energy after canonical')
    print(mf.mo_energy)

    return mf

def get_fock(mf, ref):
    '''
    mf: scf object

    return:
    scf energy
    Fock matrix

    get Fock matrix in MO basis
    set diag elements of Fock matrix as orb energy
    then set diag as zero
    this function will change the mf object
    '''

    mol = mf.mol
    
    dm = mf.make_rdm1(mf.mo_coeff, mf.mo_occ)

    h1e = mf.get_hcore() # h_uv = <u|h^core|v>
    vhf = mf.get_veff(mol, dm) # v_uv = <ui||vi>
    fock = h1e+vhf # f_uv

    e_scf = mf.energy_tot(dm, h1e, vhf)
    print('SCF energy:',e_scf)

    if ref == 'RHF':
        fock_in_mo = np.einsum('uv,up,vq->pq', fock, mf.mo_coeff, mf.mo_coeff)
        norb = len(mf.mo_energy)
        #print('difference form old orb energy and diag Fock elements:')
        #print([fock_in_mo[i,i]-mf.mo_energy[i] for i in range(0,norb)])
        print('orb energy:')
        for i in range(0,norb):
            mf.mo_energy[i] = fock_in_mo[i,i]
            fock_in_mo[i,i] = 0.0
            print(i, mf.mo_energy[i])
    elif ref == 'UHF':
        fock_in_mo = [None]*2
        fock_in_mo[0] = np.einsum('uv,up,vq->pq', fock[0], mf.mo_coeff[0], mf.mo_coeff[0])
        fock_in_mo[1] = np.einsum('uv,up,vq->pq', fock[1], mf.mo_coeff[1], mf.mo_coeff[1])
        norb = len(mf.mo_energy[0])
        #print('difference form old orb energy and diag Fock elements:')
        #print([fock_in_mo[0][i,i]-mf.mo_energy[0][i] for i in range(0,norb)])
        #print([fock_in_mo[1][i,i]-mf.mo_energy[1][i] for i in range(0,norb)])
        print('orb energy:')
        for i in range(0,norb):
            mf.mo_energy[0][i] = fock_in_mo[0][i,i]
            fock_in_mo[0][i,i] = 0.0
            mf.mo_energy[1][i] = fock_in_mo[1][i,i]
            fock_in_mo[1][i,i] = 0.0
            print('a', i, mf.mo_energy[0][i])
            print('b', i, mf.mo_energy[1][i])
    else:
        print('ERROR: reference in get_fock error',ref)
        exit()

    return e_scf, fock_in_mo

class orbital:
    '''
    mf: scf object
    mol: gto.M object
    fz: number of frozen orbital
    ref: reference state, UHF or RHF
    compGroup: computatioal group
    irrep_table: irreducible representation table
    nirrep: number of irrep
    symR: symmetry of Rx, y, and z
    dirprd: direct product table
    e_scf: scf energy
    fock: fock matrix
    mo_energy: orbital energy
    mo_occ: orbital occupation number
    mo_coeff: orbital coefficient
    mo_sym: orbital symmetry
    norb: total number of orbital after froze orbital
    norb_no_fz: tital number of orbital
    pop: population of occupied orbital for each irrep
    vir: population of vritrual orbital for each irrep
    popa: population of occupied orbital for spin a
    popb: population of occupied orbital for spin b
    vira: population of virtual orbital for spin a
    virb: population of virtual orbital for spin b
    n_occ: total number of occupied orbital
    n_vir: total number of vritual orbital

    All orbital information here are exclusive the forzed orbital are 
    resort as (irrep,energy)
    '''
    def __init__(self, mf, mol, fz=0, ref='UHF',test=False,tmpdir='./'):
        self.mf = mf # scf object
        self.mol = mol # gto.M()
        self.fz = fz # number of frozen orbital, for UHF, frozen the same number of alpha and beta orbital
        self.ref = ref # reference state, UHF or RHF
        self.test = test
        self.tmpdir = tmpdir

    def get_sym_info(self):
        '''
        get the group name, irrep, dirprd and symRx,y,z
        store in orbital class
        '''

        mol = self.mol

        # determine computatioal group: self.compGroup
        if not mol.symmetry:
            self.compGroup = 'C1'
        elif mol.groupname not in ['C1','Ci','Cs','C2','C2h','C2v','D2','D2h']:
            print('ERROR: please use Abel group for calculation')
            exit()
        else:
            self.compGroup = mol.groupname
        print('computation group:', self.compGroup)

        # determine self.irrep_table
        self.irrep_table = symm.param.IRREP_ID_TABLE[self.compGroup] 
        print('irrep_table:',self.irrep_table)

        # determine number of irrep: self.nirrep
        self.nirrep = len(self.irrep_table)
        print('number of irrep:', self.nirrep)

        # determine symmerty Rx,y,z: self.symR
        self.symR = [0, 0, 0] # C1, Ci
        if self.compGroup in ['Cs','C2','C2h']:
            self.symR[0] = 1
            self.symR[1] = 1
            self.symR[2] = 0
        elif self.compGroup in ['C2v', 'D2', 'D2h']:
            self.symR[0] = 3
            self.symR[1] = 2
            self.symR[2] = 1
        self.symR = [0,0,0]
        print('ATTENTATION: with out SOC, so we set symR=0')
        print('symRx,y,z:', self.symR)

        # determine symmerty x,y,z: self.symX
        self.symX = [0, 0, 0] # C1
        if self.compGroup == 'Ci':
            self.symX = [1, 1, 1]
        elif self.compGroup == 'Cs':
            self.symX = [0, 0, 1]
        elif self.compGroup == 'C2':
            self.symX = [1, 1, 0]
        elif self.compGroup == 'C2h':
            self.symX = [3, 3, 2]
        elif self.compGroup == 'C2v':
            self.symX = [2, 3, 0]
        elif self.compGroup == 'D2':
            self.symX = [3, 2, 1]
        elif self.compGroup == 'D2h':
            self.symX = [7, 6, 5]
        print('symX,Y,Z:', self.symX)

        self.dirprd = dirprd
        print('dirprd', self.dirprd[0:self.nirrep,0:self.nirrep])

    def resort_orb(self):
        '''
        get orbital energy, occupation, coefficient, symmerty
        population of occupied and virtual orb after resort
        store in orbital class
        sort by symmerty and energy

        only use for self.ref = UHF and RHF,
        ROHF are convert to UHF, CASSCF are convert to UHF or RHF determine by the spin
        '''
        mol = self.mol
        mf = self.mf
        fz = self.fz

        e_scf, fock = get_fock(mf, self.ref)
        self.e_scf = e_scf

        # get orbital information include mo_energy, mo_occ, mo_coeff, mo_sym, and total number norb
        # mo_coeff store as c_up
        if self.ref == 'RHF':
            self.fock = fock[fz:, fz:]
            self.mo_energy = mf.mo_energy[fz:]
            self.mo_occ = mf.mo_occ[fz:]
            self.mo_coeff = mf.mo_coeff[:, fz:]
            self.norb = len(self.mo_occ)
            self.norb_no_fz = len(mf.mo_occ)
            if not mol.symmetry or self.compGroup == 'C1':
                self.mo_sym = np.array([0]*self.norb_no_fz)
            else:
                self.mo_sym = mf.get_orbsym(mf.mo_coeff, mf.get_ovlp())[fz:]
        elif self.ref == 'UHF':
            self.fock = [0, 0]
            self.mo_energy = [0, 0]
            self.mo_occ = [0, 0]
            self.mo_coeff = [0, 0]
            self.mo_sym = [0, 0]
            for spin in [0, 1]:
                self.fock[spin] = fock[spin][fz:, fz:]
                self.mo_energy[spin] = mf.mo_energy[spin][fz:]
                self.mo_occ[spin] = mf.mo_occ[spin][fz:]
                self.mo_coeff[spin] = mf.mo_coeff[spin][:, fz:]
                self.norb = len(self.mo_occ[0])
                self.norb_no_fz = len(mf.mo_occ[0])
                if not mol.symmetry or self.compGroup == 'C1':
                    self.mo_sym[spin] = np.array([0]*self.norb_no_fz)
                else:
                    self.mo_sym[spin] = mf.get_orbsym(mf.mo_coeff, mf.get_ovlp())[spin][fz:]
        print('number of orbital:', self.norb_no_fz)
        print('number of orbital after fz:', self.norb)

        # get index of orbitals for symmetry
        if self.ref == 'RHF':
            occOrb = [i for i, ne in enumerate(self.mo_occ) if ne == 2]
            virOrb = [i for i, ne in enumerate(self.mo_occ) if ne == 0]
            # index of orbital for different irrep
            occOrbSym = []
            virOrbSym = []
            pop = [0 for _ in range(0, 8)]
            vir = [0 for _ in range(0, 8)]
            for irrep in range(0, self.nirrep):
                occOrbSym.append(
                    [i for i in occOrb if self.mo_sym[i] == irrep])
                virOrbSym.append(
                    [i for i in virOrb if self.mo_sym[i] == irrep])
                pop[irrep] = len(occOrbSym[irrep])
                vir[irrep] = len(virOrbSym[irrep])
        elif self.ref == 'UHF':
            occOrb = [0, 0]
            virOrb = [0, 0]
            for spin in [0, 1]:
                occOrb[spin] = [i for i, ne in enumerate(
                    self.mo_occ[spin]) if ne == 1]
                virOrb[spin] = [i for i, ne in enumerate(
                    self.mo_occ[spin]) if ne == 0]
            # index of orbital for different irrep
            occOrbSym = [[], []]
            virOrbSym = [[], []]
            pop = [[0 for _ in range(0, 8)], [0 for _ in range(0, 8)]]
            vir = [[0 for _ in range(0, 8)], [0 for _ in range(0, 8)]]
            for irrep in range(0, self.nirrep):
                for spin in [0, 1]:
                    occOrbSym[spin].append(
                        [i for i in occOrb[spin] if self.mo_sym[spin][i] == irrep])
                    virOrbSym[spin].append(
                        [i for i in virOrb[spin] if self.mo_sym[spin][i] == irrep])
                    pop[spin][irrep] = len(occOrbSym[spin][irrep])
                    vir[spin][irrep] = len(virOrbSym[spin][irrep])

        self.pop = pop  # population of orbital
        self.vir = vir
        print('pop', self.pop)
        print('vir', self.vir)
        if self.ref == 'RHF':
            self.popa = self.pop
            self.popb = self.pop
            self.vira = self.vir
            self.virb = self.vir
        elif self.ref == 'UHF':
            self.popa = self.pop[0]
            self.popb = self.pop[1]
            self.vira = self.vir[0]
            self.virb = self.vir[1]

        # resort orbital as the order of symmerty
        if self.ref == 'RHF':
            self.n_occ = sum(self.pop) # number of occupied orbital
            self.n_vir = sum(self.vir) # number of virtual orbital
            if self.n_occ+self.n_vir != self.norb:
                print('ERROR: number of orbital error:',self.n_occ, self.n_vir, self.norb)
            occOrb = [i for irrep in occOrbSym for i in irrep]
            virOrb = [i for irrep in virOrbSym for i in irrep]
            orb = occOrb+virOrb
            self.mo_energy = self.mo_energy[orb]
            self.mo_occ = self.mo_occ[orb]
            self.mo_sym = self.mo_sym[orb]
            self.mo_coeff = self.mo_coeff[:, orb]
            self.fock = self.fock[:, orb][orb, :]
        elif self.ref == 'UHF':
            orb = [0, 0]
            self.n_occ = [0,0]
            self.n_vir = [0,0]
            for spin in [0, 1]:
                self.n_occ[spin] = sum(self.pop[spin]) # number of occupied orbital
                self.n_vir[spin] = sum(self.vir[spin]) # number of virtual orbital
                if self.n_occ[spin]+self.n_vir[spin] != self.norb:
                    print('ERROR: number of orbital error:',self.n_occ[spin], self.n_vir[spin], self.norb)
                # change 2D list to 1D
                occOrb[spin] = [i for irrep in occOrbSym[spin] for i in irrep]
                virOrb[spin] = [i for irrep in virOrbSym[spin] for i in irrep]
                orb[spin] = occOrb[spin]+virOrb[spin]
                # resort orbital
                self.mo_energy[spin] = self.mo_energy[spin][orb[spin]]
                self.mo_occ[spin] = self.mo_occ[spin][orb[spin]]
                self.mo_sym[spin] = self.mo_sym[spin][orb[spin]]
                self.mo_coeff[spin] = self.mo_coeff[spin][:, orb[spin]]
                self.fock[spin] = self.fock[spin][:, orb[spin]][orb[spin], :]
        print('number of occupied orbital:',self.n_occ)
        print('number of virtual orbital:',self.n_vir)

        # resort by energy

        print('after sort:')
        print('-----------------------------------')
        print('energy:', self.mo_energy)
        print('occ:', self.mo_occ)
        print('symmetry:', self.mo_sym)
        print('-----------------------------------')
