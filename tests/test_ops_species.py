import bpy
from batoms import Batoms
from batoms.bio.bio import read
import pytest

def test_species():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    bpy.context.view_layer.objects.active = ch4.obj
    assert len(ch4.species) == 2
    bpy.ops.batoms.species_update()
    bpy.ops.batoms.species_add(species="O")
    assert len(ch4.species) == 3
    bpy.ops.batoms.species_remove(species="O")
    assert len(ch4.species) == 2



if __name__ == "__main__":
    test_species()
    print("\n Operator species: All pass! \n")
