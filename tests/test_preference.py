import bpy
import pytest
import logging
import importlib
import pkgutil


package = "batoms"
addon = bpy.context.preferences.addons[package]
preferences = addon.preferences


def test_enable_disable_plugin():
    from bpy.types import Collection, Object
    assert preferences.magres==True
    assert hasattr(Collection, 'Bmagres')
    preferences.magres=False
    assert not hasattr(Collection, 'Bmagres')
    preferences.magres=True
    assert hasattr(Collection, 'Bmagres')


def create_molecule():
    bpy.ops.batoms.delete()
    from batoms import Batoms
    from ase.build import molecule

    batoms = Batoms(from_ase=molecule("CO"))
    bpy.ops.batoms.delete()
    return


def read_log_lines(logger):
    """Ugly patch for reading the logger contents,
    since pytest caplog not working properly in blender (yet)
    """
    fh = logger.handlers[0]
    filename = fh.baseFilename
    with open(filename, "r") as fd:
        lines = fd.readlines()
    return lines


def import_submodules(package, recursive=True):
    """Import all submodules of a module, recursively, including subpackages
    https://stackoverflow.com/questions/3365740/how-to-import-all-submodules
    :param package: package (name or actual module)
    :type package: str | module
    :rtype: dict[str, types.ModuleType]
    """
    if isinstance(package, str):
        package = importlib.import_module(package)
    results = {}
    for loader, name, is_pkg in pkgutil.walk_packages(package.__path__):
        full_name = package.__name__ + "." + name
        results[full_name] = importlib.import_module(full_name)
        if recursive and is_pkg:
            results.update(import_submodules(full_name))
    return results


def test_child_loggers():
    """Make sure all submodules of batoms follow the setLevel rule
    https://stackoverflow.com/questions/3365740/how-to-import-all-submodules
    """
    all_submodules = import_submodules("batoms")
    preferences.logging_level = "WARNING"
    for name, mod in all_submodules.items():
        try:
            logger = mod.logger
            print(f"Test {name}.logger")
            assert logger.name == name
            assert logger.getEffectiveLevel() == logging.WARNING
        except AttributeError:
            print(f"{name}.logger not found. skip")
            continue


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
        assert logger.getEffectiveLevel() == getattr(logging, level)
    # bpy.context.preferences.addons[package].preferences.logging_level = "INFO"
    # assert(logger.level == 20)
    # bpy.context.preferences.addons[package].preferences.logging_level = "WARNING"
    # assert(logger.level == 30)


def test_logging_level_emit():
    """Test if setting logging level hierachically works
    Set the logging level to INFO, adding Batoms shows the timing info
    Set the level to WARNING and higher, timing info is supressed
    """
    root_logger = logging.getLogger("batoms")
    # INFO level reveals the timing info
    preferences.logging_level = "INFO"
    create_molecule()
    lines1 = read_log_lines(root_logger)
    assert any(["Add object" in line for line in lines1])
    # WARNING level suppresses timing info
    preferences.logging_level = "WARNING"
    create_molecule()
    lines2 = read_log_lines(root_logger)[len(lines1) :]
    assert all(["Add object" not in line for line in lines2])


if __name__ == "__main__":
    test_logging_level()
    test_logging_level_emit()
    test_child_loggers()
    print("\n Batoms: All pass! \n")
