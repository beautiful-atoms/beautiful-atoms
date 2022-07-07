import bpy
from batoms import Batoms
import pytest

try:
    from openbabel import pybel
    no_pybel = False
except ImportError:
    no_pybel = True


@pytest.mark.skipif(
    no_pybel,
    reason="Requires pybel module",
)
def test_pybel_smiles():
    bpy.ops.batoms.delete()
    smiles = 'CCO'
    mol = pybel.readstring("smi", smiles)
    mol.make3D(forcefield='mmff94', steps=100)
    mol = Batoms('c2H6O', from_pybel=mol)
    assert len(mol) == 9


@pytest.mark.skipif(
    no_pybel,
    reason="Requires pybel module",
)
def test_pybel_bond_order():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add(label = 'c6h6', formula = 'C6H6')
    c6h6 = Batoms('c6h6')
    c6h6.model_style = 1
    c6h6.bond.bond_order_auto_set()
    assert c6h6.bond.arrays['order'][0] == 2

if __name__ == "__main__":
    test_pybel_smiles()
    print("\n Pybel: All pass! \n")
