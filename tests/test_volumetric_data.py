import bpy
import pytest
from batoms.batoms import Batoms
from batoms.bio.bio import read
import numpy as np
from time import time

try:
    from _common_helpers import has_display, set_cycles_res

    use_cycles = not has_display()
except ImportError:
    use_cycles = False

extras = dict(engine="cycles") if use_cycles else {}


def test_settings():
    """key search"""
    from batoms.batoms import Batoms
    from ase.io.cube import read_cube_data
    bpy.ops.batoms.delete()
    volume, atoms = read_cube_data("../tests/datas/h2o-homo.cube")
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
    from batoms.batoms import Batoms
    from ase.io.cube import read_cube_data
    bpy.ops.batoms.delete()
    data, atoms = read_cube_data("../tests/datas/h2o-homo.cube")
    h2o = Batoms("h2o", from_ase=atoms)
    assert h2o.coll.batoms.ui_list_index_volumetric_data == 0
    h2o.volumetric_data["homo"] = data
    assert h2o.coll.batoms.ui_list_index_volumetric_data == 0
    #
    hartree, atoms = read_cube_data("../tests/datas/h2o-hartree.cube")
    h2o.volumetric_data["hartree"] = -hartree
    assert h2o.coll.batoms.ui_list_index_volumetric_data == 1
    #
    h2o.volumetric_data["diff"] = h2o.volumetric_data["hartree"] - \
        h2o.volumetric_data["homo"]
    assert h2o.coll.batoms.ui_list_index_volumetric_data == 2


def test_ops():
    """gui panel"""
    from batoms.batoms import Batoms
    from batoms.batoms import Batoms
    from ase.io.cube import read_cube_data
    bpy.ops.batoms.delete()
    data, atoms = read_cube_data("../tests/datas/h2o-homo.cube")
    h2o = Batoms("h2o", from_ase=atoms)
    assert len(h2o.volumetric_data) == 0
    bpy.ops.batoms.volumetric_data_add(
        name="hartree", filepath="../tests/datas/h2o-hartree.cube")
    assert h2o.volumetric_data.find("hartree") is not None
