import pytest


@pytest.fixture
def h2o():
    import bpy
    from batoms import Batoms
    bpy.ops.batoms.delete()
    h2o = Batoms(
        "h2o",
        species=["O", "H", "H"],
        positions=[[0, 0, 0.40], [0, -0.76, -0.2], [0, 0.76, -0.2]],
    )
    return h2o

@pytest.fixture
def c2h6so():
    import bpy
    from ase.build import molecule
    from batoms import Batoms
    bpy.ops.batoms.delete()
    c2h6so = molecule("C2H6SO")
    c2h6so = Batoms("c2h6so", from_ase=c2h6so)
    return c2h6so

@pytest.fixture
def au():
    import bpy
    from batoms import Batoms
    bpy.ops.batoms.delete()
    a = 4.08
    positions = [[0, 0, 0], [a / 2, a / 2, 0],
                 [a / 2, 0, a / 2], [0, a / 2, a / 2]]
    au = Batoms(
        label="au",
        species=["Au"] * len(positions),
        positions=positions,
        pbc=True,
        cell=(a, a, a),
    )
    return au

@pytest.fixture
def tio2():
    import bpy
    from batoms.bio.bio import read
    bpy.ops.batoms.delete()
    tio2 = read("../tests/datas/tio2.cif")
    return tio2
