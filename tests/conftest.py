import pytest


try:
    from openbabel import pybel
    has_openbabel = True
except ImportError:
    has_openbabel = False

@pytest.fixture
def h2o():
    import bpy
    from batoms import Batoms
    from ase.build import molecule
    # setup fixture h2o
    h2o = Batoms("h2o", from_ase=molecule("H2O"))
    yield h2o
    # teardown fixture h2o
    bpy.ops.batoms.delete()

@pytest.fixture
def c2h6so():
    import bpy
    from ase.build import molecule
    from batoms import Batoms
    c2h6so = molecule("C2H6SO")
    c2h6so = Batoms("c2h6so", from_ase=c2h6so)
    yield c2h6so
    bpy.ops.batoms.delete()

@pytest.fixture
def ch4():
    import bpy
    from ase.build import molecule
    from batoms import Batoms
    ch4 = molecule("CH4")
    ch4 = Batoms("ch4", from_ase=ch4)
    yield ch4
    bpy.ops.batoms.delete()

@pytest.fixture
def au():
    import bpy
    from batoms import Batoms
    from ase.build import bulk
    au = Batoms("au", from_ase=bulk("Au", cubic=True))
    yield au
    bpy.ops.batoms.delete()


@pytest.fixture
def tio2():
    import bpy
    from batoms.bio.bio import read
    tio2 = read("../tests/datas/tio2.cif")
    yield tio2
    bpy.ops.batoms.delete()
