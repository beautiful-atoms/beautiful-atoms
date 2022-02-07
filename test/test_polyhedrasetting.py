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
ch4 = molecule('CH4')
mh4 = ch4.copy()
mh4.translate([2, 2, 0])
mh4[0].symbol = 'N'
ch4 = ch4 + mh4
ch4 = Batoms('ch4', from_ase = ch4)
ch4.model_style = 2
ch4.bonds.setting[('C', 'H')].polyhedra = True
ch4.model_style = 1

def test_polyhedra():
from batoms.bio import read
from batoms.butils import removeAll
removeAll()
tio2 = read('test/datas/tio2.cif')
tio2.model_style = 2
tio2 = tio2*[3, 3, 3]
tio2.pbc = False
tio2.model_style = 2


from batoms.bio import read
from batoms.butils import removeAll
removeAll()
pk = read('test/datas/perovskite.cif')
pk.model_style = 2
pk = pk*[3, 3, 3]
pk.pbc = False
# pk.model_style = 2


def test_polyhedra():
from ase.build import molecule
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
ch4 = Batoms('ch4', from_ase = molecule('CH4'))
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