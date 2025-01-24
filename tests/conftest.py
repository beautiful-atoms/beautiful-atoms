import pytest
import bpy
import os

path = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def preferences():
    from batoms.utils.butils import get_preferences_addon

    return get_preferences_addon().preferences


@pytest.fixture
def h2o():
    from batoms import Batoms
    from ase.build import molecule

    # setup fixture h2o
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    yield h2o
    # teardown fixture h2o
    bpy.ops.batoms.delete()


@pytest.fixture
def c2h6so():
    from ase.build import molecule
    from batoms import Batoms

    c2h6so = molecule("C2H6SO")
    c2h6so = Batoms("c2h6so", from_ase=c2h6so)
    yield c2h6so
    bpy.ops.batoms.delete()


@pytest.fixture
def ch4():
    from ase.build import molecule
    from batoms import Batoms

    ch4 = molecule("CH4")
    ch4 = Batoms("ch4", from_ase=ch4)
    yield ch4
    bpy.ops.batoms.delete()


@pytest.fixture
def au():
    from batoms import Batoms
    from ase.build import bulk

    au = Batoms("au", from_ase=bulk("Au", cubic=True))
    yield au
    bpy.ops.batoms.delete()


@pytest.fixture
def tio2():
    from batoms.bio.bio import read

    tio2 = read(os.path.join(path, "datas/tio2.cif"))
    yield tio2
    bpy.ops.batoms.delete()


@pytest.fixture
def h2o_homo():
    from batoms.bio.bio import read

    h2o_homo = read(os.path.join(path, "datas/h2o-homo.cube"))
    yield h2o_homo
    bpy.ops.batoms.delete()
