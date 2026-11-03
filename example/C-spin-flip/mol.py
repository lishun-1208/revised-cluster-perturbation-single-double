from pyscf import gto, scf
from pyscf2cc import py2cc

mol = gto.M(
    atom = '''
    C 0.0  0.0  0.0''',
    basis = 'aug-cc-pVQZ',
    symmetry = 'D2h',
    charge = 0,
    spin = 2,
    verbose = 5,
)

mf = scf.UHF(mol)
mf.kernel()
mf.analyze()

py2cc.py2cc(mf, mol, 'UHF')