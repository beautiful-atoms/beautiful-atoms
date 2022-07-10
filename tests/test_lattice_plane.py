import bpy
import pytest
from ase.build import bulk
from batoms import Batoms
from batoms.bio.bio import read
import numpy as np
import os
try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}

skip_test = bool(os.environ.get("NOTEST_CUBE", 0))


def test_lattice_plane():
    bpy.ops.batoms.delete()
    au = bulk("Au", cubic=True)
    au = Batoms("au", from_ase=au)
    au.lattice_plane.settings[(1, 1, 1)] = {"distance": 3, "scale": 1.1}
    au.lattice_plane.settings[(1, 1, 0)] = {"distance": 3, "color": [0.8, 0.1, 0, 0.8]}
    print(au.lattice_plane.settings)
    au.lattice_plane.draw()
    # au.draw_cell()
    if use_cycles:
        set_cycles_res(au)
    au.get_image([0, 0, 1], output="plane-lattice-plane.png", **extras)


def test_boundary():
    # skip due to cubefile
    if skip_test:
        pytest.skip("Skip tests on cube files since $NOTEST_CUBE provided.")
    bpy.ops.batoms.delete()
    h2o = read("../tests/datas/h2o-homo.cube")
    h2o.lattice_plane.settings[(0, 0, 1)] = {"distance": 6, "boundary": True}
    h2o.lattice_plane.settings[(0, 0, -1)] = {"distance": -5, "boundary": True}
    h2o.lattice_plane.draw()
    if use_cycles:
        set_cycles_res(h2o)
    h2o.get_image([1, 0, 0], output="plane-boundary.png", **extras)


def test_lattice_plane_ops():
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(formula="Au")
    au = Batoms('Au')
    bpy.context.view_layer.objects.active = au.obj
    bpy.ops.plane.lattice_plane_add(indices=(1, 1, 1))
    au.lattice_plane.settings["1-1-1"].distance = 3.0
    bpy.ops.plane.lattice_plane_add(indices=(1, 0, 0))
    assert len(au.lattice_plane.settings) == 2
    bpy.ops.plane.lattice_plane_draw()
    bpy.ops.plane.lattice_plane_remove(name="1-0-0")
    assert len(au.lattice_plane.settings) == 1
    print(au.lattice_plane.settings)
    bpy.ops.plane.lattice_plane_draw()

def test_gui():
    """latticeplane panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(formula="Au", cubic=True)
    au = Batoms('Au')
    au.obj.select_set(True)
    assert bpy.context.scene.Blatticeplane.show == True
    bpy.context.scene.Blatticeplane.show = False
    assert au.lattice_plane.show == False

def test_gui_uilist():
    """latticeplane panel"""
    from batoms.batoms import Batoms
    bpy.ops.batoms.delete()
    bpy.ops.batoms.bulk_add(formula="Au", cubic=True)
    au = Batoms('Au')
    au.obj.select_set(True)
    assert au.coll.Blatticeplane.ui_list_index==0
    bpy.ops.plane.lattice_plane_add(indices=(1, 1, 1))
    bpy.ops.plane.lattice_plane_add(indices=(1, 0, 0))
    assert au.coll.Blatticeplane.ui_list_index==1


if __name__ == "__main__":
    test_lattice_plane()
    test_boundary()
    print("\n Lattice plane: All pass! \n")
