import pytest
from batoms.butils import removeAll
from batoms.batoms import Batoms
from batoms.bio import read
import numpy as np



def test_boundary():
    from batoms.butils import removeAll
    from batoms.batoms import Batoms
    removeAll()
    a = 3.96
    positions = [[0, 0, 0], [a/2, a/2, 0], [a/2, 0, a/2], [0, a/2, a/2]]
    pt = Batoms('pt', {'Pt': positions}, pbc = True, cell = (a, a, a))
    pt.boundary = 0.01
    atoms = pt.get_atoms_with_boundary()
    assert len(atoms) == 14
    pt.draw_cell()
    pt.get_image(output='boundary_pt.png')


def test_bond_search():
    from batoms.butils import removeAll
    from batoms.batoms import Batoms
    from batoms.bio import read
    removeAll()
    tio2 = read('datas/tio2.cif')
    tio2.boundary = 0.01
    tio2.model_style = 2
    atoms = tio2.get_atoms_with_boundary()
    assert len(atoms) == 55
    tio2.draw_cell()
    tio2.get_image(output='boundary_search.png')




if __name__ == '__main__':
    test_boundary()
    test_bond_search()
    print('\n Bondsetting: All pass! \n')