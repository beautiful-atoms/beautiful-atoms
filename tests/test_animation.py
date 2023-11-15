import bpy
from ase.io import read
from batoms import Batoms


def test_animation_molecule():
    bpy.ops.batoms.delete()
    deca = read("../tests/datas/deca_ala_md.xyz", index=":")
    deca = Batoms("deca", from_ase=deca, load_trajectory=True)
    deca.model_style = 1
    deca.nframe > 1


def test_animation_crystal():
    bpy.ops.batoms.delete()
    atoms = read("../tests/datas/tio2_10.xyz", index=":")
    tio2 = Batoms("tio2", from_ase=atoms, load_trajectory=True)
    tio2.model_style = 2
    tio2.boundary = 0.01
    tio2.nframe == 10
