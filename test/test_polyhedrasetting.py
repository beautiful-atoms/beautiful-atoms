import pytest
from batoms.butils import removeAll
from batoms.batoms import Batoms
from batoms.bio import read
import numpy as np



def test_polyhedra():
from ase.build import molecule
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
ch4 = Batoms('ch4', aseAtoms = molecule('CH4'))
ch4.bondsetting[('C', 'H')].polyhedra = True
ch4.model_style = 2
ch4.pbc = True
ch4.cell = [3, 3, 3]
ch4 = ch4*[2, 1, 1]
sel1 = ch4.selects.add('sel1', [0, 1, 2, 3, 4])
sel1.model_style = 1
    ch4.polyhedrasetting.remove('C')
    assert len(ch4.polyhedrasetting) == 1
    ch4.polyhedrasetting.add('C')
    assert len(ch4.polyhedrasetting) == 2
    ch4.polyhedrasetting['C'].color = [0.8, 0.1, 0.3, 0.3]
    ch4.model_style = 2
    ch4.render.engine = 'workbench'
    ch4.get_image([1, 1, 0], output = 'polyhedras.png')


if __name__ == '__main__':
    test_polyhedra()
    print('\n Bondsetting: All pass! \n')