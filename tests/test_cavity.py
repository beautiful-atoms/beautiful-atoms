import bpy
import pytest
from batoms.batoms import Batoms
from ase.build import molecule, bulk
from batoms.bio.bio import read
from time import time

try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}


def test_cavity():
    from batoms.bio.bio import read
    bpy.ops.batoms.delete()
    tio2 = read("../tests/datas/tio2.cif")
    tio2.boundary = 0.01
    tio2.cavity.resolution = 1
    # tio2.cavity.build_cavity()
    # tio2 *= [2, 2, 2]
    tio2.cavity.draw()
    if use_cycles:
        set_cycles_res(tio2)
    tio2.get_image([0, 1, 0], output="mof-5.png", **extras)


def test_cavity_zsm():
    from batoms.bio.bio import read
    bpy.ops.batoms.delete()
    mof = read("../tests/datas/zsm-5.cif")
    mof.boundary = 0.01
    mof.cavity.resolution = 1
    # mof.cavity.build_cavity()
    mof *= [2, 1, 1]
    mof.cavity.draw()
    if use_cycles:
        set_cycles_res(mof)
    mof.get_image([0, 1, 0], output="mof-5.png", **extras)

def test_cavity_mof():
    from batoms.bio.bio import read
    bpy.ops.batoms.delete()
    mof = read("../tests/datas/mof-5.cif")
    mof.boundary = 0.01
    mof.cavity.resolution = 1
    # mof.cavity.build_cavity()
    mof *= [2, 1, 1]
    mof.cavity.draw()
    if use_cycles:
        set_cycles_res(mof)
    mof.get_image([0, 1, 0], output="mof-5.png", **extras)

def test_cavity_ops():
    bpy.ops.batoms.delete()
    mof = read("../tests/datas/mof-5.cif")
    mof.cavity.resolution = 0.5
    bpy.context.view_layer.objects.active = mof.obj
    bpy.ops.surface.cavity_draw()
    assert len(mof.cavity.settings) == 1
    bpy.ops.surface.cavity_remove(name="0")
    print(mof.cavity.settings)
    assert len(mof.cavity.settings) == 0


if __name__ == "__main__":
    test_cavity()
    test_cavity_zsm()
    test_cavity_mof()
