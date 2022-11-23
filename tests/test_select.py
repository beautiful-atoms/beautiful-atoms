import bpy
import pytest
from ase.build import molecule, fcc111
from batoms.batoms import Batoms
import numpy as np
from time import time


def test_settings():
    """species panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    assert ch4.coll.batoms.ui_list_index_select==0
    # add
    ch4.selects.add('s1', [1])
    assert ch4.coll.batoms.ui_list_index_select==1
    # remove
    ch4.selects.remove('s1')
    assert ch4.coll.batoms.ui_list_index_select==0


def test_select():
    from ase.build import molecule, fcc111
    from batoms.batoms import Batoms
    import numpy as np
    bpy.ops.batoms.delete()
    au111 = fcc111("Au", (4, 4, 4), vacuum=0)
    au111 = Batoms("au111", from_ase=au111)
    mol = Batoms("mol", from_ase=molecule("CH3CH2OH"))
    mol.translate([5, 5, 10])
    au111 = au111 + mol
    au111.cell[2, 2] += 10
    assert len(au111.selects) == 3
    # model_style
    au111.selects["mol"].model_style = 1
    assert len(au111.bond) == 8
    assert np.isclose(au111[-1].scale, 0.4)
    au111.selects["mol"].model_style = 0
    assert np.isclose(au111[-1].scale, 1)
    assert len(au111.bond) == 0


def test_select_protein():
    from batoms.pdbparser import read_pdb
    bpy.ops.batoms.delete()
    atoms = read_pdb("../tests/datas/1ema.pdb")  # 1tim
    batoms = Batoms("protein", from_ase=atoms)
    batoms.ribbon.draw()
    sel1 = batoms.selects.add("sel1", np.where(batoms.arrays["types"] == "HETATM")[0])
    sel1.show = 1
    sel1.model_style = 1


if __name__ == "__main__":
    test_select()
    test_select_protein()
    print("\n Select: All pass! \n")
