import bpy
from batoms import Batoms
import numpy as np
from batoms.utils.butils import removeAll
import pytest

def test_batoms():
    """Create a molecule use GUI ASE"""
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.obj.select_set(True)
    bpy.context.scene.bapanel.model_style = {"1", "1", "1"}
    assert ch4.model_style[0] == 1


if __name__ == "__main__":
    test_batoms()
    print("\n GUI: All pass! \n")
