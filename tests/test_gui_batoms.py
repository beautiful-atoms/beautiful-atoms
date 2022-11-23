import bpy
from batoms import Batoms
import numpy as np
import pytest

def test_batoms():
    """Batoms panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.obj.select_set(True)
    # model_style
    assert bpy.context.scene.batoms.batoms.model_style.upper() == "BALL-AND-STICK"
    bpy.context.scene.batoms.batoms.model_style = "Space-filling"
    assert ch4.model_style == 0
    # radius_style
    assert bpy.context.scene.batoms.batoms.radius_style.upper() == "COVALENT"
    bpy.context.scene.batoms.batoms.radius_style = "VDW"
    assert ch4.radius_style == '1'
    # color_style
    assert bpy.context.scene.batoms.batoms.color_style.upper() == "JMOL"
    bpy.context.scene.batoms.batoms.color_style = "CPK"
    assert ch4.color_style == '2'
    # show
    assert bpy.context.scene.batoms.batoms.show == True
    bpy.context.scene.batoms.batoms.show = False
    assert ch4.show[0] == False
