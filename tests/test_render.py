import bpy
import pytest
import numpy as np
from ase.build import molecule, fcc111
from batoms.render import Render
from batoms import Batoms

try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}


def test_render():
    bpy.ops.batoms.delete()
    render = Render("test")
    assert isinstance(render, Render)
    render.studiolight = "paint.sl"
    render.viewport = [1, 0, 0]
    render.resolution = [200, 200]
    # from collection
    render = Render("test")
    assert isinstance(render, Render)
    assert render.studiolight == "paint.sl"
    assert render.resolution == [200, 200]


def test_render_gpu():
    bpy.ops.batoms.delete()
    h2o = molecule("H2O")
    h2o = Batoms("h2o", from_ase=h2o)
    h2o.render.init()
    h2o.render.gpu = True


def test_render_init():
    bpy.ops.batoms.delete()
    h2o = molecule("H2O")
    h2o = Batoms("h2o", from_ase=h2o)
    h2o.render.init()
    h2o.render.viewport = [1, 0, 0]
    h2o.render.resolution = [200, 200]
    if use_cycles:
        set_cycles_res(h2o)
    h2o.get_image(**extras)


def test_render_setter():
    bpy.ops.batoms.delete()
    h2o = molecule("H2O")
    h2o = Batoms("h2o", from_ase=h2o)
    h2o.render.init()
    h2o.render.viewport = [1, 0, 0]
    nh3 = molecule("NH3")
    nh3 = Batoms("nh3", from_ase=nh3)
    nh3.translate([3, 3, 0])
    nh3.render = Render("nh3")
    nh3.render.viewport = [0, 1, 0]
    nh3.render.lights.add("1", type="SUN", direction=[1, 0, 0])
    nh3.render.resolution = [200, 200]
    if use_cycles:
        set_cycles_res(nh3)
    nh3.get_image(output="render-setter.png", **extras)


def test_render_au111():
    bpy.ops.batoms.delete()
    atoms = fcc111("Au", size=(4, 4, 4), vacuum=0)
    au111 = Batoms(label="au111", from_ase=atoms)
    au111.cell[2, 2] += 5
    au111.render.resolution = [200, 200]
    if use_cycles:
        set_cycles_res(au111)
    au111.get_image(**extras)
    au111.render.engine = "eevee"
    au111.get_image([1, 1, 1], output="au111-eevee.png", **extras)
    au111.render.engine = "cycles"
    au111.get_image([1, 1, 1], output="au111-cycles.png", **extras)


if __name__ == "__main__":
    test_render()
    test_render_setter()
    test_render_au111()
    print("\n render: All pass! \n")
