import bpy

# assume on github test
try:
    from _common_helpers import use_cycles

    update = not use_cycles()
except ImportError:
    update = False


def test_use_batoms_startup():
    bpy.ops.batoms.delete()
    # no display, assume on github test
    if update:
        bpy.ops.batoms.use_batoms_startup()


def test_update():
    bpy.ops.batoms.delete()
    # no display, assume on github test
    if update:
        bpy.ops.batoms.update()


def test_pip_install_package():
    bpy.ops.batoms.delete()
    # no display, assume on github test
    if update:
        bpy.ops.batoms.pip_install_package()


def test_use_batoms_preference():
    bpy.ops.batoms.delete()
    # no display, assume on github test
    if update:
        bpy.ops.batoms.use_batoms_preference()
