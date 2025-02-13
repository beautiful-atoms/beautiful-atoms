import bpy
import numpy as np
import os

path = os.path.dirname(os.path.abspath(__file__))

try:
    from _common_helpers import use_cycles

    use_cycles = not use_cycles()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}


def test_settings():
    """key search"""
    from batoms.batoms import Batoms
    from ase.io.cube import read_cube_data

    bpy.ops.batoms.delete()
    volume, atoms = read_cube_data(os.path.join(path, "datas/h2o-homo.cube"))
    h2o = Batoms("h2o", from_ase=atoms, volume={"homo": volume})
    assert h2o.volumetric_data.bpy_setting["homo"].shape[:] == volume.shape
    assert h2o.volumetric_data["homo"].shape == volume.shape
    assert len(h2o.volumetric_data) == 1
    h2o.volumetric_data["electrostatic"] = volume + 1
    assert len(h2o.volumetric_data) == 2
    assert h2o.volumetric_data.find("homo") is not None
    assert h2o.volumetric_data.find("electrostatic") is not None
    h2o.volumetric_data.remove("homo")
    assert h2o.volumetric_data.find("homo") is None


def test_gui():
    """gui panel"""
    from batoms.batoms import Batoms
    from ase.io.cube import read_cube_data

    bpy.ops.batoms.delete()
    data, atoms = read_cube_data(os.path.join(path, "datas/h2o-homo.cube"))
    h2o = Batoms("h2o", from_ase=atoms)
    assert h2o.coll.batoms.ui_list_index_volumetric_data == 0
    h2o.volumetric_data["homo"] = data
    assert h2o.coll.batoms.ui_list_index_volumetric_data == 0
    #
    hartree, atoms = read_cube_data(os.path.join(path, "datas/h2o-hartree.cube"))
    h2o.volumetric_data["hartree"] = -hartree
    assert h2o.coll.batoms.ui_list_index_volumetric_data == 1
    #
    h2o.volumetric_data["diff"] = (
        h2o.volumetric_data["hartree"] - h2o.volumetric_data["homo"]
    )
    assert h2o.coll.batoms.ui_list_index_volumetric_data == 2


def test_ops_add():
    """gui panel"""
    from batoms.batoms import Batoms
    from ase.io.cube import read_cube_data

    bpy.ops.batoms.delete()
    data, atoms = read_cube_data(os.path.join(path, "datas/h2o-homo.cube"))
    h2o = Batoms("h2o", from_ase=atoms)
    assert len(h2o.volumetric_data) == 0
    bpy.context.view_layer.objects.active = h2o.obj
    bpy.ops.batoms.volumetric_data_add(
        name="hartree", filepath=os.path.join(path, "datas/h2o-hartree.cube")
    )
    assert h2o.volumetric_data.find("hartree") is not None
    bpy.ops.batoms.volumetric_data_remove(name="hartree")
    assert h2o.volumetric_data.find("hartree") is None


def test_ops_create(h2o_homo):
    from ase.io.cube import read_cube_data

    h2o = h2o_homo
    bpy.context.view_layer.objects.active = h2o.obj
    hartree, atoms = read_cube_data(os.path.join(path, "datas/h2o-hartree.cube"))
    h2o.volumetric_data["hartree"] = -hartree
    bpy.ops.batoms.volumetric_data_create(
        name="diff", select_data1="h2o_homo", select_data2="h2o_homo", operator="Minus"
    )
    assert h2o.coll.batoms.ui_list_index_volumetric_data == 2
    assert np.isclose(h2o.volumetric_data["diff"], 0).all()


def test_repeat():
    from batoms.bio.bio import read

    bpy.ops.batoms.delete()
    h2o = read(os.path.join(path, "datas/h2o-homo.cube"))
    shape = h2o.volumetric_data["h2o_homo"].shape
    M = [2, 2, 1]
    h2o._volumetric_data *= M
    new_shape = h2o.volumetric_data["h2o_homo"].shape
    assert list(new_shape) == [shape[i] * M[i] for i in range(3)]
