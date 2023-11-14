import bpy
import pytest
from ase.build import bulk
from batoms import Batoms
from ase.build import graphene_nanoribbon


@pytest.mark.skip(reason="Need to run manualy by mouse")
def test_force_field():
    bpy.ops.batoms.delete()
    gnr = graphene_nanoribbon(2, 2, type="armchair", saturated=True, vacuum=3.5)
    gnr.pbc = False
    gnr = Batoms("gnr", from_ase=gnr)


@pytest.mark.skip(reason="Need to run manualy by mouse")
def test_force_field_al():
    bpy.ops.batoms.delete()
    al = Batoms("al", from_ase=bulk("Al"))
    al = al * [1, 20, 1]
    al.pbc = False
