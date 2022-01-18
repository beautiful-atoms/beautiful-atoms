import pytest
from batoms.butils import removeAll
from batoms.batoms import Batoms
from batoms.bio import read
import numpy as np



def test_bondsetting():
from ase.io import read
from ase.build import molecule
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
mof = read('test/datas/mof-5.cif')
h2o = molecule('CH3CH2OH')
h2o.pbc = True
h2o.center(3)
h2o = Batoms('h2o', from_ase = mof)
# h2o.bondsetting.remove('C-C')
# h2o.bondsetting.remove('O-O')
h2o.model_style = 1
bpy.data.objects['h2o_bond_center'].data.attributes['order'].data[0].value = 2
bpy.data.objects['h2o_bond_center'].data.attributes['style'].data[0].value = 2

def test_model_style():

from batoms.batoms import Batoms
from ase.build import molecule
from batoms.butils import removeAll
removeAll()
mol = molecule('H2O')
h2o = Batoms('h2o', from_ase = mol)
h2o.pbc = True
h2o.cell = [3, 3, 3]
h2o = h2o*[2, 1, 1]
indices = [0, 1, 2]
sel1 = h2o.selects.add('sel1', indices)
sel0 = h2o.selects['sel0']
sel1.model_style = 1
sel0.scale = 0.5
sel1.scale = [2, 2, 2]
h2o.draw_space_filling()
    h2o.replace('H', 'H_1', [0])
    h2o.bondsetting.remove(['O', 'H_1'])
    h2o['H_1'].color = [0.1, 0.8, 0, 1.0]
    h2o.bondsetting.add(['O', 'H_1'])
    h2o.model_style = 1
    h2o.render.viewport=[1, 0, 0]
    h2o.render.engine = 'eevee'
    h2o.get_image(output = 'bonds-replace.png')

def test_add():
    from batoms.build import bulk
    from batoms.butils import removeAll
    removeAll()
    au = bulk('au', 'Au')
    au = au*[2, 2, 2]
    assert len(au.bondsetting) == 0
    au.bondsetting.add(['Au', 'Au'])
    assert len(au.bondsetting) == 1
    print('Pass add!')

def test_polyhedra():
    from batoms.bio import read
    from batoms.butils import removeAll
    removeAll()
    tio2 = read('datas/tio2.cif')
    tio2.boundary = 0.01
    tio2.bondsetting.remove(('Ti', 'O'))
    assert len(tio2.bondsetting) == 1
    tio2.bondsetting.add(('Ti', 'O'))
    assert len(tio2.bondsetting) == 2

def test_search_bond():
    from batoms.bio import read
    from batoms.butils import removeAll
    removeAll()
    pk = read('datas/perovskite.cif')
    pk.repeat([2, 2, 2])
    pk.boundary = 0.01
    pk.model_style = 2
    pk.draw_cell()
    pk.get_image([0, 1, 0], engine = 'eevee', output = 'perovskite.png')

def test_search_bond_2():
    from batoms.bio import read
    from batoms.butils import removeAll
    removeAll()
    mol = read('datas/anthraquinone.cif')
    mol.boundary = 0.01
    mol.draw_cell()
    mol.model_style = 1
    mol.get_image([1, -0.3, 0.1], engine = 'eevee', output = 'anthraquinone.png')

def test_search_bond_urea():
    from batoms.bio import read
    from batoms.butils import removeAll
    removeAll()
    mol = read('datas/urea.cif')
    mol.boundary = 2
    mol.draw_cell()
    mol.model_style = 1
    mol.get_image([1, -0.3, 0.1], engine = 'eevee', output = 'urea.png')

def test_search_bond_3():
    """
    Big system
    """
    from batoms.bio import read
    from batoms.butils import removeAll
    removeAll()
    mof = read('datas/mof-5.cif')
    mof.boundary = 0.01
    mof.bondsetting[('Zn', 'O')].polyhedra = True
    mof.model_style = 1
    mof.get_image([0, 1, 0], engine = 'eevee', output = 'mof-5.png')

def test_high_order_bond():
    """
    High order bond
    """
    from ase.build import molecule
    from batoms import Batoms
    from batoms.butils import removeAll
    removeAll()
    c6h6 = molecule('C6H6')
    c6h6 = Batoms('c6h6', atoms = c6h6)
    for i in range(6):
        c6h6.replace('C', 'C_%s'%i, [0])

    c6h6.bondsetting[('C_1', 'C_0')].order = 2
    c6h6.bondsetting[('C_3', 'C_2')].order = 2
    c6h6.bondsetting[('C_5', 'C_4')].order = 2
    c6h6.model_style = 1
    c6h6.get_image([0, 0, 1], engine = 'eevee', output = 'c6h6.png')

def test_hydrogen_bond():
    from ase.build import molecule
    from batoms.batoms import Batoms
    from batoms.butils import removeAll
    removeAll()
    h2o = molecule('H2O')
    h2o2 = molecule('H2O')
    h2o2.rotate(90, 'x')
    h2o2.translate([0, 0, 3])
    h2o = h2o + h2o2
    h2o = Batoms(label = 'h2o', atoms = h2o)
    h2o.bondsetting[('H', 'O')].min = 2.0
    h2o.bondsetting[('H', 'O')].max = 3.0
    h2o.bondsetting[('H', 'O')].width = 0.01
    h2o.bondsetting[('H', 'O')].style = '2'
    h2o.model_style = 1
    h2o.get_image([1, 0 ,0], engine = 'eevee', output = 'bonds-hb.png')

if __name__ == '__main__':
    test_replace()
    test_add()
    test_polyhedra()
    test_hydrogen_bond()
    test_search_bond()
    test_search_bond_2()
    test_search_bond_3()
    test_high_order_bond()
    print('\n Bondsetting: All pass! \n')