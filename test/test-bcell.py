import pytest
from batoms.cell import Bcell
import numpy as np

def test_bcell():
    """
    """
    cell = Bcell(label = 'pt', array = [2, 2, 2])
    assert isinstance(cell, Bcell)
    # properties
    cell[2, 2] += 5
    assert cell[2, 2] == 7
    # repeat
    cell.repeat([2, 2, 2])
    assert np.allclose(cell.array, np.array([[4, 0, 0], [0, 4, 0], [0, 0, 14]]))
    cell2 = cell.copy('pt-2')
    assert isinstance(cell2, Bcell)

if __name__ == '__main__':
    test_bcell()
    print('\n Bcell: All pass! \n')