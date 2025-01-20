"""This module is only to be used together with tests
"""
import os


def use_cycles():
    """Use environment variable to set Cycles as a boolean.

    export USE_CYCLES=1  # Enable Cycles
    export USE_CYCLES=0  # Disable Cycles

    github action
    env:
        USE_CYCLES: "1"
    """
    return os.environ.get("USE_CYCLES", "0").lower() in ("1", "true", "yes")


# Very bad resolution only for cycles output
CYCLES_TEST_RES = [20, 20]


def set_cycles_res(obj, res=CYCLES_TEST_RES):
    """Set the cycles render resolution, default to 20x20"""
    assert hasattr(obj, "render")
    obj.render.resolution = res
    return
