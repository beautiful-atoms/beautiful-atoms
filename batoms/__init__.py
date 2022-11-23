bl_info = {
    "name": "Batoms toolbar",
    "author": "Xing Wang",
    "version": (2, 2, 0),
    "blender": (3, 0, 0),
    "location": "File -> Import -> Batoms (xyz, cif, pdb, ...)",
    "description": """Python module for drawing and
rendering atomic structures using blender.""",
    "warning": "",
    "category": "Import-Export",
    "doc_url": "https://beautiful-atoms.readthedocs.io/en/latest/",
    "tracker_url": "https://github.com/beautiful-atoms/beautiful-atoms/issues/new/choose",
}

from time import time
tstart0 = time()
import bpy
# install pip dependencies
from .install import pip_dependencies
pip_dependencies.install()

from batoms.batoms import Batoms



from . import (
    logger,
    preferences,
    plugins,
    internal_data,
    pip_dependencies,
    ops,
    gui,
    console,
)

logger.set_logger(bl_info["version"])



import importlib

modules = ['bond',
            'polyhedra',
            'render',
            'ribbon']

def enable_module():
    for key in modules:
        module = importlib.import_module("batoms.{}".format(key))
        module.register_class()

def disable_module():
    for key in modules:
        module = importlib.import_module("batoms.{}".format(key))
        module.unregister_class()


def register():
    # dependencies
    pip_dependencies.register_class()
    preferences.register_class()
    # class
    internal_data.register_class()
    # class
    ops.register_class()
    gui.register_class()
    # manual
    ops.register_manual_map()
    # menu
    ops.register_menu()
    gui.register_menu()
    # keymap
    gui.register_keymap()
    # hook
    console.register_hook()
    # modules
    enable_module()
    # plugins
    plugins.enable_plugin()
    logger.root_logger.info("Batoms init time: {:.2f}".format(time() - tstart0))
    logger.update_logging_level()



def unregister():
    # dependencies
    pip_dependencies.unregister_class()
    # class
    internal_data.unregister_class()
    ops.unregister_class()
    gui.unregister_class()
    # manual
    ops.unregister_manual_map()
    # menu
    ops.unregister_menu()
    gui.unregister_menu()
    # keymap
    gui.unregister_keymap()
    # hook
    console.unregister_hook()
    disable_module()
    plugins.disable_plugin()
    preferences.unregister_class()

if __name__ == "__main__":

    register()
