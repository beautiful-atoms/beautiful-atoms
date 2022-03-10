import bpy
from batoms import Batoms
import numpy as np
from batoms.utils.butils import removeAll
import pytest

def test_ase_molecule():
    """Create a molecule use GUI ASE"""
    removeAll()
    bpy.context.scene.asepanel.formula = "NH3"
    bpy.context.scene.asepanel.label = "nh3"
    bpy.ops.batoms.add_molecule()
    nh3 = Batoms('nh3')
    assert len(nh3) == 4

def test_ase_bulk():
    """Create an bulk use GUI ASE"""
    removeAll()
    bpy.context.scene.asepanel.formula = "Au"
    bpy.context.scene.asepanel.label = "au"
    bpy.ops.batoms.add_bulk()
    au = Batoms('au')
    assert len(au) == 1

def test_ase_surface():
    """Create an surface use GUI ASE"""
    removeAll()
    bpy.context.scene.asepanel.formula = "Au"
    bpy.context.scene.asepanel.label = "au"
    bpy.ops.batoms.add_surface()
    au = Batoms('au')
    assert len(au) == 1


if __name__ == "__main__":
    test_ase_molecule()
    test_ase_bulk()
    print("\n Batoms: All pass! \n")
