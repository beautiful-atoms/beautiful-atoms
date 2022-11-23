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

def test_edit_bonds():
    """Edit bonds panel"""
    from batoms.batoms import Batoms
    import numpy as np
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.bond.obj.data.vertices.foreach_set('select', [1, 0, 0, 0, 0])
    ch4.bond[0].order = 2
    # order
    ch4.obj.select_set(False)
    bpy.context.view_layer.objects.active = ch4.bond.obj
    assert np.isclose(bpy.context.scene.Bbond.order, ch4.bond[0].order)
    ch4.bond.obj.data.vertices.foreach_set('select', [1, 0, 0, 0])
    bpy.context.scene.Bbond.order = 1
    assert np.isclose(ch4.bond[0].order, bpy.context.scene.Bbond.order)
