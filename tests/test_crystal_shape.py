import bpy
import pytest
from ase.build import bulk
from batoms import Batoms
import numpy as np
import os
try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}
skip_test = bool(os.environ.get("NOTEST_CUBE", 0))


def test_crystal_shape():
    from batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(label = 'au', formula = "Au", cubic=True)
    au = Batoms("au")
    au.crystal_shape.settings[(1, 1, 1)] = {
        "distance": 3, "crystal": True, "symmetry": True}
    au.crystal_shape.settings[(0, 0, 1)] = {'distance': 3,
    'crystal': True, 'symmetry': True}
    print(au.crystal_shape.settings)
    au.crystal_shape.draw(origin=au.cell.center)
    if use_cycles:
        set_cycles_res(au)
    au.get_image([0, 0, 1], output="plane-crystal.png", **extras)


def test_crystal_shape_ops():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(formula="Au", cubic=True)
    au = Batoms('Au')
    bpy.context.view_layer.objects.active = au.obj
    bpy.ops.plane.crystal_shape_add(indices=(1, 1, 1))
    au.crystal_shape.settings["1-1-1"].symmetry = True
    au.crystal_shape.settings["1-1-1"].distance = 3.0
    bpy.ops.plane.crystal_shape_add(indices=(1, 0, 0))
    assert len(au.crystal_shape.settings) == 2
    bpy.ops.plane.crystal_shape_draw()
    bpy.ops.plane.crystal_shape_remove(name="1-0-0")
    assert len(au.crystal_shape.settings) == 8
    print(au.crystal_shape.settings)
    bpy.ops.plane.crystal_shape_draw()

def test_settings():
    """key search"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(formula="Au", cubic=True)
    au = Batoms('Au')
    au.obj.select_set(True)
    au.crystal_shape.settings.add((1, 1, 1))
    assert au.crystal_shape.settings.find((1, 1, 1)) is not None
    assert au.crystal_shape.settings.find('1-1-1') is not None
    au.crystal_shape.settings.remove((1, 1, 1))
    assert au.crystal_shape.settings.find((1, 1, 1)) is None

def test_gui():
    """crystal shape panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(formula="Au", cubic=True)
    au = Batoms('Au')
    au.obj.select_set(True)
    assert bpy.context.scene.Bcrystalshape.show == True
    bpy.context.scene.Bcrystalshape.show = False
    assert au.crystal_shape.show == False

def test_crystalshape_uilist():
    """crystalshape panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(formula="Au", cubic=True)
    au = Batoms('Au')
    au.obj.select_set(True)
    assert au.crystal_shape.settings.ui_list_index==0
    bpy.ops.plane.crystal_shape_add(indices=(1, 1, 1))
    bpy.ops.plane.crystal_shape_add(indices=(1, 0, 0))
    assert au.coll.Bcrystalshape.ui_list_index==1


if __name__ == "__main__":
    test_crystal_shape()
    print("\n Crystal shape: All pass! \n")
