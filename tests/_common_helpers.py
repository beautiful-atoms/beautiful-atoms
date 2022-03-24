"""This module is only to be used together with tests
"""
import os


def has_display():
    """Detect display, should work in linux systems; windows not tested"""
    display_string = os.environ.get("DISPLAY", "")
    return len(display_string) > 0


# Very bad resolution only for cycles output
CYCLES_TEST_RES = [20, 20]


def set_cycles_res(obj, res=CYCLES_TEST_RES):
    """Set the cycles render resolution, default to 20x20"""
    assert hasattr(obj, "render")
    obj.render.resolution = res
    return
