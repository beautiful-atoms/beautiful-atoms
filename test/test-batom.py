import pytest
from batoms.batom import Batom
import numpy as np
from batoms.butils import removeAll

def test_batom():
    """
    """
    removeAll()
    positions = np.array([[0, 0, 0], [1.52, 0, 0]])
    h = Batom('h2o', 'H', positions = positions)
    assert isinstance(h, Batom)
    h.translate([0, 0, 2])
    assert np.allclose(h.positions, np.array([[0, 0, 2], [1.52, 0, 2]]))
    # properties
    h.scale = 2
    assert np.allclose(h.scale, np.array([2, 2, 2]))
    #
    h1 = Batom('1 o', 'H_1', [[0, 0, 2], [2, 0, 2]])
    h = h + h1
    h2 = h.copy('h2o', 'H_2')
    #
    assert isinstance(h2, Batom)
    h3=Batom('atom_h2o_H_2')
    assert isinstance(h3, Batom)

def test_positions():
    from batoms.batom import Batom
    import numpy as np
    removeAll()
    positions = np.array([[0, 0, 0], [1.52, 0, 0]])
    h = Batom('h2o', 'H', positions = positions)
    npositions = positions - np.array([0, 0, 5])
    h.positions = npositions


def test_batom_animation():
    from batoms.batom import Batom
    import numpy as np
    removeAll()
    positions = np.array([[0, 0 ,0], [1.52, 0, 0]])
    o = Batom('co2', 'O', positions)
    images = []
    for i in range(10):
        images.append(positions + np.array([i, 0, 0]))
    o.load_frames(images)


if __name__ == '__main__':
    test_batom()
    test_positions()
    test_batom_animation()
    print('\n Batom: All pass! \n')