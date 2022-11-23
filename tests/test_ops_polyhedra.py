import bpy
from batoms import Batoms
from batoms.bio.bio import read
import pytest

def test_polyhedra():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    bpy.context.view_layer.objects.active = ch4.obj
    assert len(ch4.polyhedra.settings) == 2
    bpy.ops.batoms.polyhedra_draw()
    bpy.ops.batoms.polyhedra_remove(species="C")
    assert len(ch4.polyhedra.settings) == 1
    bpy.ops.batoms.polyhedra_add(species="C")
    assert len(ch4.polyhedra.settings) == 2


if __name__ == "__main__":
    test_polyhedra()
    print("\n Operator polyhedra: All pass! \n")
