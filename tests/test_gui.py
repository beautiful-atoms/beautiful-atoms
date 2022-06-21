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

def test_cell():
    """Cell panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.cell = [1, 2, 3]
    ch4.obj.select_set(True)
    # model_style
    assert np.isclose(bpy.context.scene.batoms.cell.cell_a0, 1)
    bpy.context.scene.batoms.cell.cell_a0 = 3
    assert np.isclose(ch4.cell[0, 0], 3)

def test_render():
    """Render panel"""
    from batoms.batoms import Batoms
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
    assert np.isclose(bpy.context.scene.batoms.render.viewport_x, 1)
    bpy.context.scene.batoms.render.viewport_y = 2
    assert ch4.render.viewport[1] == 2
    # camera
    ch4.render.camera.type = "ORTHO"
    assert bpy.context.scene.batoms.render.camera_type == "ORTHO"
    # Some GUI tests can only be tested manually by opening Blender
    # bpy.context.scene.batoms.render.camera_type = "PERSP"
    # assert ch4.render.camera.type == 'PERSP'
    ch4.render.camera.ortho_scale = 20
    assert np.isclose(bpy.context.scene.batoms.render.ortho_scale, 20)
    # bpy.context.scene.batoms.render.ortho_scale = 5
    # assert np.isclose(ch4.render.camera.ortho_scale, 5)
    # light
    ch4.render.lights["Default"].direction = [1, 0, 0]
    assert np.isclose(bpy.context.scene.batoms.render.light_direction_x, 1)
    ch4.render.lights["Default"].energy = 20
    assert np.isclose(bpy.context.scene.batoms.render.energy, 20)
    # bpy.context.scene.batoms.render.energy = 5
    # assert np.isclose(ch4.render.lights["Default"].energy, 5)

def test_batom():
    """Batom panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.obj.select_set(True)
    # scale
    assert np.isclose(bpy.context.scene.batoms.batom.scale, ch4[0].scale)
    bpy.context.scene.batoms.batom.scale = 1
    assert np.isclose(ch4[0].scale, bpy.context.scene.batoms.batom.scale)

def test_bond():
    """Bond panel"""
    from batoms.batoms import Batoms
    import numpy as np
    bpy.ops.batoms.delete()
    bpy.ops.batoms.molecule_add()
    ch4 = Batoms('CH4')
    ch4.bonds.obj.data.vertices.foreach_set('select', [1, 0, 0, 0, 0])
    ch4.bonds[0].order = 2
    # order
    ch4.obj.select_set(False)
    bpy.context.view_layer.objects.active = ch4.bonds.obj
    assert np.isclose(bpy.context.scene.batoms.bond.order, ch4.bonds[0].order)
    ch4.bonds.obj.data.vertices.foreach_set('select', [1, 0, 0, 0])
    bpy.context.scene.batoms.bond.order = 1
    assert np.isclose(ch4.bonds[0].order, bpy.context.scene.batoms.bond.order)


if __name__ == "__main__":
    test_batoms()
    test_cell()
    test_render()
    test_batom()
    test_bond()
    print("\n GUI: All pass! \n")
