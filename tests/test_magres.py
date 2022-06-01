import bpy
import pytest
from batoms.batoms import Batoms
from ase.io import read
import numpy as np
from time import time


def test_ms():
    from batoms import Batoms
    from ase.io import read
    atoms = read('../tests/datas/ethanol.magres')
    ms_array = atoms.get_array('ms')
    ethanol = Batoms("ethanol", from_ase=atoms)
    ethanol.model_style = 1
    ethanol.magres.tensors = ms_array
    ethanol.magres.setting['1'].scale = 0.005
    ethanol.magres.draw()


if __name__ == "__main__":
    test_ms()
    print("\n Magres: All pass! \n")
