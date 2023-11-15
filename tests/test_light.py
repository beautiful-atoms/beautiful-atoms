import bpy
from batoms.render import Render, Lights
from batoms.build import molecule

try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}


def test_lights():
    """ """
    bpy.ops.batoms.delete()
    r = Render(label="h2o")
    assert isinstance(r.lights, Lights)
    # properties
    r.lights.add("1")
    r.lights["1"].type = "POINT"
    assert r.lights["1"].type == "POINT"
    r.lights.add("2")
    r.lights.remove("2")
    assert len(r.lights) == 2


def test_light_direction():
    bpy.ops.batoms.delete()
    h2o = molecule("h2o", "H2O")
    h2o.render.lights["Default"].type = "POINT"
    h2o.render.lights["Default"].energy = 1000
    assert h2o.render.lights["Default"].type == "POINT"
    if use_cycles:
        set_cycles_res(h2o)
    h2o.get_image([0, 0, 1], output="light-direction.png", **extras)
