from pyscf import gto, scf
from pyscf2cc import py2cc

mol = gto.M(
    atom = '''
    O 0.0  0.00000000 -0.13209669
    H 0.0  1.43152878  0.97970006
    H 0.0 -1.43152878  0.97970006''',
    unit='B',
    basis = 'aug-cc-pVTZ',
    symmetry = 'C2v',
    charge = 0,
    spin = 0,
    verbose = 5,
)

mf = scf.RHF(mol)
mf.kernel()
mf.analyze()

py2cc.py2cc(mf, mol, 'RHF')