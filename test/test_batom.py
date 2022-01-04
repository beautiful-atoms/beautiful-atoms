import pytest
from batoms.batom import Batom
import numpy as np
from batoms.butils import removeAll

def test_batom():
    """
    """
    from batoms import Batom
    import numpy as np
    from batoms.butils import removeAll
    removeAll()
    positions = np.array([[0, 0, 0], [1.52, 0, 0]])
    h = Batom('atom_H', species = 'H', positions = positions)
    assert isinstance(h, Batom)
    #
    # from object
    h2=Batom('atom_H')
    assert isinstance(h2, Batom)
    # copy
    h3 = h.copy('atom_H3')
    assert isinstance(h3, Batom)
    #
    h.translate([0, 0, 2])
    assert np.allclose(h.positions, np.array([[0, 0, 2], [1.52, 0, 2]]))
    # properties
    h.radius_style = 'VDW'
    print('radius: ', h.radius)
    assert np.allclose(h.radius, 1.20)
    h.scale = 2
    assert np.allclose(h.scale, np.array([2, 2, 2]))
    assert np.allclose(h.size, 2.40)
    # extend
    h1 = Batom('atom_H_1', species = 'H_1', positions = [[0, 0, 2], [2, 0, 2]])
    h4 = h + h1
    assert isinstance(h4, Batom)


def test_positions():
    from batoms import Batom
    import numpy as np
    from batoms.butils import removeAll
    removeAll()
    positions = np.array([[0, 0, 0], [1.52, 0, 0]])
    h = Batom('H', positions = positions)
    npositions = positions - np.array([0, 0, 5])
    h.positions = npositions
    assert np.allclose(h.positions, np.array([[0, 0, -5], [1.52, 0, -5]]))

def test_batom_animation():
    from batoms import Batom
    import numpy as np
    from batoms.butils import removeAll
    removeAll()
    positions = np.array([[0, 0 ,0], [1.52, 0, 0]])
    o = Batom('O', positions = positions)
    images = []
    for i in range(10):
        images.append(positions + np.array([i, 0, 0]))     
    
    o.frames = images
    # repeat
    cell = np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]])
    o.repeat([2, 1, 1], cell)
    # join
    positions = np.array([[0.76, 0 ,0]])
    o_1 = Batom('O_1', positions = positions)
    images = []
    for i in range(10):
        images.append(positions + np.array([i, 0, 0]))
    
    o_1.frames = images
    o.extend(o_1)
    # delete
    o.delete([4])

def test_occupancy():
    """
    """
    from batoms import Batom
    import numpy as np
    from batoms.butils import removeAll
    removeAll()
    positions = np.array([[0, 0, 0], [3, 0, 0]])
    o = Batom('o', elements = {'O': 0.833, 'N': 0.167}, positions = positions)
    assert len(o.elements) == 2


if __name__ == '__main__':
    test_batom()
    test_positions()
    test_batom_animation()
    test_occupancy()
    print('\n Batom: All pass! \n')