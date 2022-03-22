import bpy
import pytest
from batoms import Batoms
import numpy as np
from batoms.render import Camera
from ase.build import fcc111

try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}


def test_camera():
    """ """
    bpy.ops.batoms.delete()
    camera = Camera(label="h2o", name="1")
    assert isinstance(camera, Camera)
    # properties
    camera.type = "PERSP"
    assert camera.type == "PERSP"
    camera.lens = 50
    assert camera.lens == 50
    camera.target = [5, 0, 0]


def test_camera_ortho():
    bpy.ops.batoms.delete()
    au111 = fcc111("Au", (4, 4, 4), vacuum=0)
    au111 = Batoms("au111", from_ase=au111)
    au111.cell[2, 2] += 10
    # au111.draw_cell()
    if use_cycles:
        set_cycles_res(au111)
    au111.get_image([1, 0, 0], output="camera_ortho.png", **extras)
    au111.get_image(
        [0, 0, 1], center=au111.positions[-1], output="camera-center.png", **extras
    )


def test_camera_persp():
    bpy.ops.batoms.delete()
    au111 = fcc111("Au", (4, 10, 4), vacuum=0)
    au111 = Batoms("au111", from_ase=au111)
    au111.cell[2, 2] += 10
    # au111.draw_cell()
    au111.render.camera.type = "PERSP"
    if use_cycles:
        set_cycles_res(au111)
    au111.get_image([1, 0, 0], output="camera-persp.png", **extras)


if __name__ == "__main__":
    test_camera()
    test_camera_ortho()
    test_camera_persp()
    print("\n camera: All pass! \n")
