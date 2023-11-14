def test_slicebonds(ch4):
    """Setting sliced Bonds"""
    import numpy as np

    ch4.model_style = 1
    # one atoms
    assert ch4.bond[0].order == 1
    ch4.bond[1].order = 2
    assert ch4.bond[1].order == 2
    ch4.bond[1].show = 0
    assert ch4.bond[1].show == 0
    # > one bonds
    assert np.isclose(ch4.bond[0:2].order, np.array([1, 2])).all()
    ch4.bond[0:2].order = 2
    assert np.isclose(ch4.bond[0:2].order, np.array([2, 2])).all()
