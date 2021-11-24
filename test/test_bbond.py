import pytest
from batoms.batom import Batom
import numpy as np
from batoms.butils import removeAll

def test_bbond():
    """
    """
    from batoms.bond import Bbond
    import numpy as np
    from batoms.butils import removeAll

    removeAll()
    positions = np.array([[0, 0, 0, 0, 0, 1, 1], [1.52, 0, 0, 0, 0, 1, 1]])
    h = Bbond('h2o', 'O_H_H', positions=positions, 
                battr_inputs = {'bbond':{'flag': True, 'segments': 16}})
    assert isinstance(h, Bbond)
    h.translate([0, 0, 2])
    assert len(h) == 2


def test_bbond_animation():
    """
    """
    from batoms.bond import Bbond
    import numpy as np
    from batoms.butils import removeAll
    removeAll()
    positions = np.array([[0, 0, 0, 0, 0, 1, 1], [1.52, 0, 0, 0, 0, 1, 1]])
    h = Bbond('h2o', 'O_H_H', positions=positions,
            battr_inputs = {'bbond':{'flag': True, 'segments': 16}})
    images = []
    for i in range(10):
        images.append(positions + np.array([i, 0, 0, 0, 0, 0, i*0.1]))

    h.frames = images


if __name__ == '__main__':
    test_bbond()
    test_bbond_animation()
    print('\n Batom: All pass! \n')