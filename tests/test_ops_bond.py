import bpy
from batoms import Batoms
from batoms.bio.bio import read
import pytest

def test_bond():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    bpy.context.view_layer.objects.active = ch4.obj
    bpy.ops.bond.bond_pair_remove(name="C-H")
    assert len(ch4.bonds.setting) == 1
    bpy.ops.bond.bond_pair_add(species1="C", species2="H")
    assert len(ch4.bonds.setting) == 2
    bpy.ops.bond.draw()

if __name__ == "__main__":
    test_bond()
    print("\n Operator bond: All pass! \n")