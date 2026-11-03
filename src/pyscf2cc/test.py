import numpy as np
from pyscf import gto, scf, mcscf
import symm_info
import py2cc 

mol = gto.M(
    atom='''
    O 0. 0. 0.
    O 0. 0. 1.''',
    symmetry='D2h',
    basis='cc-pVDZ',
    charge=1,
    spin=1,
    verbose=5
)

mf = scf.UHF(mol)
mf.kernel()
mf.analyze()

#casscf = mcscf.CASSCF(mf, 4, 4).run()
#print('------------------')
#print(casscf.mol.spin)
#print(casscf.mo_energy)
#print(mf.mo_energy)

#mo_occ = [np.zeros(60),np.zeros(60)]
#mo_occ[0][0:8]=1
#mo_occ[1][0:7]=1
#mo_occ[1][8]=1
print('------ start pyscf2cc ------')
orb,path=py2cc.py2cc(mf, mol, 'UHF', fz=0, canonical=False)
py2cc.get_2e_int_outcore(orb,path)
