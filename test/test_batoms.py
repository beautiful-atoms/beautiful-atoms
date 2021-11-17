import pytest
from batoms.batom import Batom
from batoms.batoms import Batoms
from batoms.butils import removeAll
import numpy as np


def test_batom():
    removeAll()
    h = Batom(label = 'h2o', species = 'H', positions = [[0, -0.76, -0.2], [0, 0.76, -0.2]])
    o = Batom(label = 'h2o', species = 'O', positions = [[0, 0, 0.40]])
    h2o = Batoms('h2o', [h, o])
    assert isinstance(h2o, Batoms)

def test_batoms():
    """
    """
    import math
    from batoms.butils import removeAll
    from batoms import Batoms
    removeAll()
    h2o = Batoms('h2o', {'O': [[0, 0, 0.40]], 'H': [[0, -0.76, -0.2], [0, 0.76, -0.2]]})
    assert isinstance(h2o, Batoms)
    # properties
    h2o.cell = [3, 3, 3]
    h2o.pbc = True
    h2o.model_type = 1
    h2o.translate([0, 0, 2])
    h2o.model_type = 1
    h2o.rotate(math.pi/2.0, 'Z')
    #
    h2o_2 = h2o.copy('h2o_2')
    assert isinstance(h2o_2, Batoms)
    h2o_3 = Batoms(label='h2o_2')
    #
    # delete
    del h2o_3['H'][[0]]
    assert len(h2o_3['H']) == 1
    #
    h2o_3.replace('O', 'C', [0])
    #
    h2o['O'][0] = [0, 0, 2.0]
    index = [0, 1]
    h2o['H'][index][:, 0] += 2

def test_from_coll():
    from batoms.butils import removeAll
    from batoms import Batoms
    removeAll()
    batoms = Batoms('h2o', {'O': [[0, 0, 0.40]], 'H': [[0, -0.76, -0.2], [0, 0.76, -0.2]]})
    h2o = Batoms('h2o')
    assert isinstance(h2o, Batoms)
    h2o.render.run([0, 0, 1], output = 'batoms_from_coll.png')
    # properties


def test_ase_species():
    from batoms.butils import removeAll
    from batoms import Batoms
    removeAll()
    from ase.build import molecule
    h2o = molecule('H2O')
    h2o.new_array('species', np.array(h2o.get_chemical_symbols(), dtype = 'U20'))
    h2o.arrays['species'][1] = 'H_1'
    h2o.arrays['species'][2] = 'H_test2'
    h2o = Batoms('h2o', atoms = h2o)
    assert h2o.species == ['H_1', 'H_test2', 'O']

def test_set_positions():
    from batoms.butils import removeAll
    from batoms import Batoms
    from ase.build import molecule
    removeAll()
    h2o = Batoms('h2o', {'O': [[0, 0, 0.40]], 'H': [[0, -0.76, -0.2], [0, 0.76, -0.2]]})
    mol = molecule('H2O')
    mol.positions[:, 2] += 5
    h2o.positions = mol

def test_extend():
    from ase.build import molecule
    from batoms.butils import removeAll
    removeAll()
    co = molecule('CO')
    co = Batoms(label = 'co', atoms = co)
    co.translate([0, 0, 2])
    h2o = molecule('H2O')
    h2o = Batoms(label = 'h2o', atoms = h2o)
    h2o.extend(co)

def test_transform():
    from batoms.butils import removeAll
    from batoms.bio import read
    removeAll()
    tio2 = read('test/datas/tio2.cif')
    tio2_t = tio2.transform([[1, 1, 0, 0], [-1, 1, 0, 0], [0, 0, 1, 0]])
    assert len(tio2_t) == 12

def test_repeat():
    from ase.build import molecule
    from batoms.butils import removeAll
    removeAll()
    h2o = molecule('H2O')
    h2o = Batoms(label = 'h2o', atoms = h2o)
    h2o.cell = [3, 3, 3]
    h2o.pbc = True
    h2o.repeat([2, 2, 2])
    h2o.model_type = 1

def test_repeat_animation():
    from batoms.butils import removeAll
    from batoms.bio import read
    removeAll()
    tio2 = read('datas/tio2_10.xyz', index = ':')
    tio2.set_frames()   
    tio2.repeat([2, 2, 2])
    tio2.model_type = 1

def test_canvas():
    from ase.build import fcc111
    from batoms.batoms import Batoms
    from batoms.butils import removeAll
    removeAll()
    atoms = fcc111('Pt', size = (4, 4, 4), vacuum=0)
    pt111 = Batoms(label = 'pt111', atoms = atoms)
    pt111.cell[2, 2] += 5
    pt111.render([1, 1, 1])

def test_boundary():
    from batoms.bio import read
    from batoms.butils import removeAll
    removeAll()
    mof = read('datas/mof-5.cif')
    mof.boundary = 0.01
    mof.boundary = 1

def test_cavity():
    from batoms.bio import read
    from batoms.butils import removeAll
    removeAll()
    mof = read('datas/mof-5.cif')
    mof.draw_cavity_sphere(9.0, boundary = [[0.2, 0.8], [0.2, 0.8], [0.2, 0.8]])
    mof.model_type = 2
    mof.draw_cell()
    mof.render.light_energy = 5
    mof.render.run([1, 0, 0], output = 'cavity.png')


def test_get_angles():
    from ase.build import molecule
    from batoms.butils import removeAll
    removeAll()
    atoms = molecule('H2O')
    h2o = Batoms(atoms = atoms, label = 'h2o')
    angle = h2o.get_angle('H', 0, 'O', 0, 'H', 1)
    d = h2o.get_distances('H', 0, 'H', 1)



if __name__ == '__main__':
    test_batom()
    test_batoms()
    test_from_coll()
    test_ase_species()
    test_set_positions()
    test_cavity()
    test_get_angles()
    test_repeat()
    test_repeat_animation()
    print('\n Batoms: All pass! \n')