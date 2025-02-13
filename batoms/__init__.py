from .batoms import Batoms  # noqa: E402

__version__ = "2.3.0"
__all__ = ["Batoms"]

from . import (  # noqa: E402
    logger,
    preferences,
    plugins,
    internal_data,
    ops,
    gui,
    console,
)

logger.set_logger(__version__)


def enable_module():
    from .bond import register_class as register_bond_class
    from .polyhedra import register_class as register_polyhedra_class
    from .render import register_class as register_render_class
    from .ribbon import register_class as register_ribbon_class

    register_bond_class()
    register_polyhedra_class()
    register_render_class()
    register_ribbon_class()


def disable_module():
    from .bond import unregister_class as unregister_bond_class
    from .polyhedra import unregister_class as unregister_polyhedra_class
    from .render import unregister_class as unregister_render_class
    from .ribbon import unregister_class as unregister_ribbon_class

    unregister_bond_class()
    unregister_polyhedra_class()
    unregister_render_class()
    unregister_ribbon_class()


def register():
    from time import time

    tstart0 = time()
    # dependencies
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
