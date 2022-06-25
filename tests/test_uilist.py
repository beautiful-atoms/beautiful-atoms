import bpy
from batoms import Batoms
from batoms.bio.bio import read
import pytest


def test_species():
    """species panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.obj.select_set(True)
    assert ch4.coll.batoms.ui_list_index_species==0
    bpy.ops.batoms.species_add(species='H_1')
    assert ch4.coll.batoms.ui_list_index_species==2



if __name__=="__main__":
    print("\n UIList: All pass! \n")
