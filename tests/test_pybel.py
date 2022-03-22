from batoms import Batoms
import pytest

try:
    from openbabel import pybel
except ImportError:
    pybel = False


@pytest.mark.skipif(
    pybel,
    reason="Requires pybel module",
)
def test_pybel_smiles():
    smiles = 'CCO'
    mol = pybel.readstring("smi", smiles)
    mol.make3D(forcefield='mmff94', steps=100)
    mol = Batoms('c2H6O', from_pybel=mol)
    assert len(mol) == 9

if __name__ == "__main__":
    test_pybel_smiles()
    print("\n Pybel: All pass! \n")
