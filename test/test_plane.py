import pytest
from ase.build import bulk
from batoms import Batoms
from batoms.utils.butils import removeAll
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
    removeAll()
    au = bulk("Au", cubic=True)
    au = Batoms("au", from_ase=au)
    au.planesetting[(1, 1, 1)] = {"distance": 3, "scale": 1.1}
    au.planesetting[(1, 1, 0)] = {"distance": 3, "color": [0.8, 0.1, 0, 0.8]}
    au.planesetting.draw_lattice_plane()
    # au.draw_cell()
    if use_cycles:
        set_cycles_res(au)
    au.get_image([0, 0, 1], output="plane-lattice-plane.png", **extras)


def test_crystal_shape():
    removeAll()
    au = bulk("Au", cubic=True)
    au = Batoms("au", from_ase=au)
    au.planesetting[(1, 1, 1)] = {"distance": 3, "crystal": True, "symmetry": True}
    # au.planesetting[(0, 0, 1)] = {'distance': 3, 'crystal': True, 'symmetry': True}
    au.planesetting.draw_crystal_shape(origin=au.cell.center)
    if use_cycles:
        set_cycles_res(au)
    au.get_image([0, 0, 1], output="plane-crystal.png", **extras)


def test_boundary():
    # skip due to cubefile
    if skip_test:
        pytest.skip("Skip tests on cube files since $NOTEST_CUBE provided.")
    removeAll()
    h2o = read("/home/xing/ase/batoms/h2o-homo.cube")
    h2o.planesetting[(0, 0, 1)] = {"distance": 6, "boundary": True}
    h2o.planesetting[(0, 0, -1)] = {"distance": -5, "boundary": True}
    h2o.planesetting.draw_lattice_plane()
    if use_cycles:
        set_cycles_res(h2o)
    h2o.get_image([1, 0, 0], output="plane-boundary.png", **extras)


if __name__ == "__main__":
    test_lattice_plane()
    test_crystal_shape()
    test_boundary()
    print("\n Bcell: All pass! \n")
