bl_info = {
    "name": "Batoms toolbar",
    "author": "Xing Wang",
    "version": (2, 1, 0),
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
    custom_property,
    pip_dependencies,
    ops,
    gui,
    console,
)

logger.set_logger(bl_info["version"])




def register():
    # dependencies
    pip_dependencies.register_class()
    preferences.register_class()
    # class
    custom_property.register_class()
    # class
    ops.register_class()
    gui.register_class()
    if bpy.context.preferences.addons['batoms'].preferences.real_module:
        from . import modal
        modal.register_class()
    # manual
    ops.register_manual_map()
    # menu
    ops.register_menu()
    gui.register_menu()
    # keymap
    gui.register_keymap()
    # hook
    console.register_hook()
    logger.root_logger.info("Batoms init time: {:.2f}".format(time() - tstart0))
    logger.update_logging_level()



def unregister():
    # dependencies
    pip_dependencies.unregister_class()
    # class
    custom_property.unregister_class()
    ops.unregister_class()
    gui.unregister_class()
    if bpy.context.preferences.addons['batoms'].preferences.real_module:
        from . import modal
        modal.unregister_class()
    # manual
    ops.unregister_manual_map()
    # menu
    ops.unregister_menu()
    gui.unregister_menu()
    # keymap
    gui.unregister_keymap()
    # hook
    console.unregister_hook()
    preferences.unregister_class()

if __name__ == "__main__":

    register()
