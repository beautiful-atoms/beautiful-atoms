import bpy
from batoms import Batoms
from batoms.bio.bio import read
import pytest

def test_bond_settings():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    bpy.context.view_layer.objects.active = ch4.obj
    bpy.ops.bond.bond_pair_remove(name="C-H")
    assert len(ch4.bond.settings) == 1
    bpy.ops.bond.bond_pair_add(species1="C", species2="H")
    assert len(ch4.bond.settings) == 2
    bpy.ops.bond.draw()

def test_bond_hydrogen_bond():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add(label = 'ch3oh', formula='CH3OH')
    ch3oh = Batoms('ch3oh')
    assert len(ch3oh.bond) == 5
    bpy.ops.bond.show_hydrogen_bond()
    assert len(ch3oh.bond) == 8

def test_bond_search_bond():
    from batoms.bio.bio import read
    bpy.ops.batoms.delete()
    tio2 = read("../tests/datas/tio2.cif")
    tio2.boundary = 0.01
    tio2.model_style = 1
    assert tio2.bond.show_search == False
    bpy.ops.bond.show_search()
    assert tio2.bond.show_search == True



if __name__ == "__main__":
    test_bond_settings()
    print("\n Operator bond: All pass! \n")
