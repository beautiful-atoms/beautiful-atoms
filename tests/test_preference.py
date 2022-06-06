import bpy
import pytest


package = "batoms"

def test_logging_level():
    """logging level
    """
    bpy.context.preferences.addons[package].preferences.logging_level = "DEBUG"
    from batoms.batoms import logger
    assert(logger.level == 10)
    bpy.context.preferences.addons[package].preferences.logging_level = "INFO"
    assert(logger.level == 20)
    bpy.context.preferences.addons[package].preferences.logging_level = "WARNING"
    assert(logger.level == 30)


if __name__ == "__main__":
    test_logging_level()
    print("\n Batoms: All pass! \n")
