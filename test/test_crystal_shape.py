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
    bpy.ops.batoms.delete()
    au = bulk("Au", cubic=True)
    au = Batoms("au", from_ase=au)
    au.crystal_shape.setting[(1, 1, 1)] = {
        "distance": 3, "crystal": True, "symmetry": True}
    au.crystal_shape.setting[(0, 0, 1)] = {'distance': 3, 
    'crystal': True, 'symmetry': True}
    au.crystal_shape.draw(origin=au.cell.center)
    if use_cycles:
        set_cycles_res(au)
    au.get_image([0, 0, 1], output="plane-crystal.png", **extras)


if __name__ == "__main__":
    test_crystal_shape()
    print("\n Crystal shape: All pass! \n")
