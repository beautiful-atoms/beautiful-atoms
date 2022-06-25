import bpy
import pytest
from batoms.batoms import Batoms
from batoms.bio.bio import read
import numpy as np
from time import time

try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}
import os

skip_test = bool(os.environ.get("NOTEST_CUBE", 0))

def test_slice():
    if skip_test:
        pytest.skip("Skip tests on cube files since $NOTEST_CUBE provided.")
    bpy.ops.batoms.delete()
    h2o = read("../tests/datas/h2o-homo.cube")
    h2o.isosurface.settings["1"] = {"level": -0.001}
    h2o.isosurface.settings["2"] = {"level": 0.001, "color": [0, 0, 0.8, 0.5]}
    h2o.isosurface.draw()
    h2o.lattice_plane.settings[(1, 0, 0)] = {"distance": 6, "slicing": True}
    h2o.lattice_plane.draw()
    if use_cycles:
        set_cycles_res(h2o)
    h2o.get_image([0, 0, 1], **extras)


def test_diff():
    if skip_test:
        pytest.skip("Skip tests on cube files since $NOTEST_CUBE provided.")
    bpy.ops.batoms.delete()
    h2o = read("../tests/datas/h2o-homo.cube", label = "h2o")
    volume = h2o.volume
    h2o.volume = volume + 0.1
    assert np.allclose(h2o.volume, volume + 0.1)
    h2o.isosurface.settings[1].level = 0.008
    h2o.isosurface.settings[2] = {"level": -0.008, "color": [0, 0, 1, 0.8]}
    h2o.model_style = 1
    h2o.isosurface.draw()
    if use_cycles:
        set_cycles_res(h2o)
    else:
        h2o.render.resolution = [200, 200]
    h2o.get_image([0, 0, 1], output="h2o-homo-diff-top.png", **extras)
    h2o.get_image([1, 0, 0], output="h2o-homo-diff-side.png", **extras)


def test_isosurface_ops():
    bpy.ops.batoms.delete()
    h2o = read("../tests/datas/h2o-homo.cube")
    bpy.context.view_layer.objects.active = h2o.obj
    bpy.ops.surface.isosurface_draw()
    assert len(h2o.isosurface.settings) == 1
    bpy.ops.surface.isosurface_remove(name="1")
    print(h2o.isosurface.settings)
    assert len(h2o.isosurface.settings) == 0
    bpy.ops.surface.isosurface_add(name="1")
    assert len(h2o.isosurface.settings) == 1
    bpy.ops.surface.isosurface_draw()


def test_isosurface_uilist():
    """isosurface panel"""
    bpy.ops.batoms.delete()
    h2o = read("../tests/datas/h2o-homo.cube")
    bpy.context.view_layer.objects.active = h2o.obj
    h2o.obj.select_set(True)
    assert h2o.coll.Bisosurface.ui_list_index==0
    bpy.ops.surface.isosurface_add(name="2")
    assert h2o.coll.Bisosurface.ui_list_index==1


if __name__ == "__main__":
    test_slice()
    test_diff()
    print("\n Isosurface.setting: All pass! \n")
