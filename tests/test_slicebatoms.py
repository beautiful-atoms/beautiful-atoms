def test_slicebatoms(ch4):
    """Setting sliced Batoms"""
    import numpy as np

    # one atoms
    ch4.model_style = 1
    assert ch4[0].species == "C"
    ch4[1].scale = 1
    assert ch4[1].scale == 1
    ch4[1].show = 0
    assert ch4[1].show == 0
    #     ch4[1].species = "Cl"
    #     assert ch4[1].species == "Cl"
    assert np.isclose(
        ch4[1].position, np.array([0.62911803, 0.62911803, 0.62911803])
    ).all()
    # > one atoms
    assert ch4[0:2].species[0] == "C"
    print(ch4[:].scale)
    assert np.isclose(ch4[0:2].scale, np.array([0.4, 1])).all()
    ch4[0:2].scale = 2
    assert np.isclose(ch4[0:2].scale, np.array([2, 2])).all()
