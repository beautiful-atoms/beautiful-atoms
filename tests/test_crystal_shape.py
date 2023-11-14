import bpy
import pytest
import os

try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}
skip_test = bool(os.environ.get("NOTEST_CUBE", 0))


def test_crystal_shape(au):
    au.crystal_shape.settings[(1, 1, 1)] = {
        "distance": 3,
        "crystal": True,
        "symmetry": True,
    }
    au.crystal_shape.settings[(0, 0, 1)] = {
        "distance": 3,
        "crystal": True,
        "symmetry": True,
    }
    au.crystal_shape.draw(origin=au.cell.center)
    if use_cycles:
        set_cycles_res(au)
    au.get_image([0, 0, 1], output="plane-crystal.png", **extras)


def test_crystal_shape_ops(au):
    # need spglib to get the symmetry
    pytest.importorskip("spglib")
    bpy.context.view_layer.objects.active = au.obj
    bpy.ops.plane.crystal_shape_add(indices=(1, 1, 1))
    au.crystal_shape.settings["1-1-1"].symmetry = True
    au.crystal_shape.settings["1-1-1"].distance = 3.0
    bpy.ops.plane.crystal_shape_add(indices=(1, 0, 0))
    assert len(au.crystal_shape.settings) == 2
    bpy.ops.plane.crystal_shape_draw()
    bpy.ops.plane.crystal_shape_remove(name="1-0-0")
    assert len(au.crystal_shape.settings) == 8
    bpy.ops.plane.crystal_shape_draw()


def test_settings(au):
    """key search"""
    bpy.context.view_layer.objects.active = au.obj
    au.crystal_shape.settings.add((1, 1, 1))
    assert au.crystal_shape.settings.find((1, 1, 1)) is not None
    assert au.crystal_shape.settings.find("1-1-1") is not None
    au.crystal_shape.settings.remove((1, 1, 1))
    assert au.crystal_shape.settings.find((1, 1, 1)) is None


def test_gui(au):
    """crystal shape panel"""
    bpy.context.view_layer.objects.active = au.obj
    assert bpy.context.scene.Bcrystalshape.show is True
    bpy.context.scene.Bcrystalshape.show = False
    assert au.crystal_shape.show is False


def test_crystalshape_uilist(au):
    """crystalshape panel"""
    bpy.context.view_layer.objects.active = au.obj
    assert au.crystal_shape.settings.ui_list_index == 0
    bpy.ops.plane.crystal_shape_add(indices=(1, 1, 1))
    bpy.ops.plane.crystal_shape_add(indices=(1, 0, 0))
    assert au.coll.Bcrystalshape.ui_list_index == 1
