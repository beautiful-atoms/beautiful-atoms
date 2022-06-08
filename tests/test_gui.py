import bpy
from batoms import Batoms
import numpy as np
import pytest

def test_batoms():
    """Batoms panel"""
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.obj.select_set(True)
    # model_style
    assert bpy.context.scene.bapanel.model_style.upper() == "BALL-AND-STICK"
    bpy.context.scene.bapanel.model_style = "Space-filling"
    assert ch4.model_style[0] == 0
    # radius_style
    assert bpy.context.scene.bapanel.radius_style.upper() == "COVALENT"
    bpy.context.scene.bapanel.radius_style = "VDW"
    assert ch4.radius_style == '1'
    # color_style
    assert bpy.context.scene.bapanel.color_style.upper() == "JMOL"
    bpy.context.scene.bapanel.color_style = "CPK"
    assert ch4.color_style == '2'
    # show
    assert bpy.context.scene.bapanel.show == True
    bpy.context.scene.bapanel.show = False
    assert ch4.show[0] == False

def test_cell():
    """Cell panel"""
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.cell = [1, 2, 3]
    ch4.obj.select_set(True)
    # model_style
    assert np.isclose(bpy.context.scene.clpanel.cell_a0, 1)
    bpy.context.scene.clpanel.cell_a0 = 3
    assert np.isclose(ch4.cell[0, 0], 3)

def test_render():
    """Render panel"""
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.cell = [1, 2, 3]
    ch4.obj.select_set(True)
    # bpy.ops.batoms.render_add()
    ch4.render.init()
    assert ch4._render is not None
    # viewport
    ch4.render.viewport = [1, 0, 0]
    assert np.isclose(bpy.context.scene.repanel.viewport_x, 1)
    bpy.context.scene.repanel.viewport_y = 2
    assert ch4.render.viewport[1] == 2
    # camera
    assert bpy.context.scene.repanel.camera_type == "ORTHO"
    bpy.context.scene.repanel.camera_type = "PERSP"
    assert ch4.render.camera.type == 'PERSP'
    ch4.render.camera.ortho_scale = 20
    assert np.isclose(bpy.context.scene.repanel.ortho_scale, 20)
    bpy.context.scene.repanel.ortho_scale = 5
    assert np.isclose(ch4.render.camera.ortho_scale, 5)
    # light
    ch4.render.lights["Default"].direction = [1, 0, 0]
    assert np.isclose(bpy.context.scene.repanel.light_direction_x, 1)
    ch4.render.lights["Default"].energy = 20
    assert np.isclose(bpy.context.scene.repanel.light_energy, 20)
    bpy.context.scene.repanel.light_energy = 5
    assert np.isclose(ch4.render.lights["Default"].energy, 5)


if __name__ == "__main__":
    test_batoms()
    test_cell()
    test_render()
    print("\n GUI: All pass! \n")
