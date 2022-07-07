import bpy
from batoms import Batoms
from batoms.bio.bio import read
import pytest


    
def test_bondpair():
    """bondpair panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.obj.select_set(True)
    assert ch4.coll.Bbond.ui_list_index==1
    bpy.ops.bond.bond_pair_add(species1='H', species2='H')
    assert ch4.coll.Bbond.ui_list_index==2
