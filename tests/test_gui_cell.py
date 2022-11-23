import bpy
from batoms import Batoms
from batoms.bio.bio import read
import pytest


def test_cell():
    """Cell panel"""
    from batoms.batoms import Batoms
    import numpy as np
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.obj.select_set(True)
    # model_style
    assert np.isclose(bpy.context.scene.batoms.cell.cell_a0, 0)
    bpy.context.scene.batoms.cell.cell_a0 = 3
    assert np.isclose(bpy.context.scene.batoms.cell.cell_a0, 3)
    assert np.isclose(ch4.cell[0, 0], 3)
