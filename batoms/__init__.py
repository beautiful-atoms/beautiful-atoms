bl_info = {
    "name": "Batoms toolbar",
    "author": "Xing Wang",
    "version": (2, 0, 0),
    "blender": (3, 0, 0),
    "location": "File -> Import -> Batoms (xyz, cif, pdb, ...)",
    "description": """Python module for drawing and
rendering atomic structures using blender.""",
    "warning": "",
    "category": "Import-Export",
    "doc_url": "https://beautiful-atoms.readthedocs.io/en/latest/",
    "tracker_url": "https://github.com/superstar54/beautiful-atoms/issues/new/choose",
}

# install pip dependencies
from . import install_dependencies
install_dependencies.install()

from batoms.batoms import Batoms


from . import (
    custom_property,
    install_dependencies,
    ops,
    gui,
    console,
    modal,
)



def register():
    # class
    custom_property.register_class()
    # class
    ops.register_class()
    gui.register_class()
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


def unregister():
    # class
    custom_property.unregister_class()
    ops.unregister_class()
    gui.unregister_class()
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

if __name__ == "__main__":

    register()
