import bpy
import pytest
import logging


package = "batoms"
addon = bpy.context.preferences.addons[package]
preferences = addon.preferences

def create_molecule():
    bpy.ops.batoms.delete()
    from batoms import Batoms
    from ase.build import molecule
    batoms = Batoms(from_ase=molecule("CO"))
    bpy.ops.batoms.delete()


def test_logging_level():
    """logging level. 
    Since now all the child loggers are "delegation to the parent" 
    https://docs.python.org/3/library/logging.html#logging.Logger.setLevel, 
    explicitly testing logger.level will not work. Instead, use getEffectiveLevel to check 
    the effect of propagation
    """
    from batoms.batoms import logger

    for level in ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]:
        preferences.logging_level = level
        assert(logger.getEffectiveLevel() == getattr(logging, level))
    # bpy.context.preferences.addons[package].preferences.logging_level = "INFO"
    # assert(logger.level == 20)
    # bpy.context.preferences.addons[package].preferences.logging_level = "WARNING"
    # assert(logger.level == 30)

def test_logging_level_emit(caplog):
    """Test if setting logging level hierachically works
    """
    # caplog.set_level(logging.DEBUG)
    create_molecule()
    root_logger = logging.getLogger("batoms")
    # root_logger.

    
    


if __name__ == "__main__":
    test_logging_level()
    print("\n Batoms: All pass! \n")
