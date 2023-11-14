import bpy


def test_magres():
    from batoms import Batoms
    from ase.io import read

    bpy.ops.batoms.delete()
    atoms = read("../tests/datas/ethanol.magres")
    ms_array = atoms.arrays.pop("ms")
    print(ms_array.shape)
    for i in range(3):
        for j in range(3):
            atoms.set_array("ms_%d%d" % (i, j), ms_array[:, i, j])

    efg_array = atoms.arrays.pop("efg")
    for i in range(3):
        for j in range(3):
            atoms.set_array("efg_%d%d" % (i, j), efg_array[:, i, j])

    ethanol = Batoms("ethanol", from_ase=atoms)
    ethanol.model_style = 1
    ethanol.magres.settings["1"].scale = 0.005
    ethanol.magres.draw()


def test_magres_uilist():
    """magres panel"""
    from batoms import Batoms
    from ase.io import read

    bpy.ops.batoms.delete()
    atoms = read("../tests/datas/ethanol.magres")
    ethanol = Batoms("ethanol", from_ase=atoms)
    bpy.context.view_layer.objects.active = ethanol.obj
    assert ethanol.coll.Bmagres.ui_list_index == 0
    bpy.ops.surface.magres_add(name="2")
    assert ethanol.coll.Bmagres.ui_list_index == 1
