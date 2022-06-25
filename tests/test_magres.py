import bpy
import pytest
from batoms.batoms import Batoms
from ase.io import read
import numpy as np
from time import time


def test_magres():
    from batoms import Batoms
    from ase.io import read
    bpy.ops.batoms.delete()
    atoms = read('../tests/datas/ethanol.magres')
    ethanol = Batoms("ethanol", from_ase=atoms)
    ethanol.model_style = 1
    ethanol.magres.settings['1'].scale = 0.005
    ethanol.magres.draw()

def test_magres_uilist():
    """magres panel"""
    from batoms import Batoms
    from ase.io import read
    bpy.ops.batoms.delete()
    atoms = read('../tests/datas/ethanol.magres')
    ethanol = Batoms("ethanol", from_ase=atoms)
    assert ethanol.coll.Bmagres.ui_list_index==0
    bpy.ops.surface.magres_add(name="2")
    assert ethanol.coll.Bmagres.ui_list_index==1


if __name__ == "__main__":
    test_magres()
    print("\n Magres: All pass! \n")
