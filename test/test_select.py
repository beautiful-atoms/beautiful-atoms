import bpy
import pytest
from ase.build import molecule, fcc111
from batoms.utils.butils import removeAll
from batoms.batoms import Batoms
import numpy as np
from time import time


def test_select():
    bpy.ops.batoms.delete()
    au111 = fcc111("Au", (4, 4, 4), vacuum=0)
    au111 = Batoms("au111", from_ase=au111)
    mol = Batoms("mol", from_ase=molecule("CH3CH2OH"))
    mol.translate([5, 5, 10])
    au111 = au111 + mol
    au111.cell[2, 2] += 10
    assert len(au111.selects) == 3
    au111.selects["mol"].model_style = 1


def test_select_protein():
    from batoms.pdbparser import read_pdb
    bpy.ops.batoms.delete()
    atoms = read_pdb("datas/1ema.pdb")  # 1tim
    batoms = Batoms("protein", from_ase=atoms)
    batoms.ribbon.draw()
    sel1 = batoms.selects.add("sel1", np.where(batoms.arrays["types"] == "HETATM")[0])
    sel1.show = 1
    sel1.model_style = 1


if __name__ == "__main__":
    test_select()
    test_select_protein()
    print("\n Select: All pass! \n")
