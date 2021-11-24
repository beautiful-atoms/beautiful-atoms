import pytest
from batoms.batom import Batom
from batoms.batoms import Batoms
from batoms.butils import removeAll
import numpy as np


def test_from_batom():
    removeAll()
    h = Batom(label = 'h2o', species = 'H', positions = [[0, -0.76, -0.2], [0, 0.76, -0.2]])
    o = Batom(label = 'h2o', species = 'O', positions = [[0, 0, 0.40]])
    h2o = Batoms('h2o', [h, o])
    assert isinstance(h2o, Batoms)
    assert len(h2o.species) == 2
    assert len(h2o) == 3

def test_batoms():
    """
    """
    import math
    from batoms.butils import removeAll
    from batoms import Batoms
    removeAll()
    h2o = Batoms('h2o', {'O': [[0, 0, 0.40]], 'H': [[0, -0.76, -0.2], [0, 0.76, -0.2]]})
    assert isinstance(h2o, Batoms)
    assert len(h2o.species) == 2
    assert len(h2o) == 3
    # properties
    h2o.cell = [3, 3, 3]
    h2o.pbc = True
    assert h2o.pbc
    h2o.model_type = 1
    assert h2o.model_type == 1
    h2o.translate([0, 0, 2])
    assert np.allclose(h2o['O'].positions, np.array([[0, 0, 2.4]]))
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
    assert len(h2o.species) == 2
    assert len(h2o) == 3
    h2o.translate([2, 0, 0])

def test_render():
    from batoms.butils import removeAll
    from batoms import Batoms
    removeAll()
    h2o = Batoms('h2o', {'O': [[0, 0, 0.40]], 'H': [[0, -0.76, -0.2], [0, 0.76, -0.2]]})
    h2o.get_image(output = 'batoms_render.png')
    h2o.get_image(viewport = [1, 0, 0], output = 'batoms_render_viewport.png')
    h2o.get_image(canvas=[5, 5, 5], output = 'batoms_render.png_canvas')

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
    species = h2o.species
    species.sort()
    assert species == ['H_1', 'H_test2', 'O']

def test_get_and_set_positions():
    from batoms.butils import removeAll
    from batoms import Batoms
    from ase.build import molecule
    removeAll()
    h2o = Batoms('h2o', {'O': [[0, 0, 0.40]], 'H': [[0, -0.76, -0.2], [0, 0.76, -0.2]]})
    positions = h2o.positions
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
    assert len(h2o.species) == 3
    assert len(h2o) == 5

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
    assert len(h2o) == 24

def test_repeat_animation():
    from batoms.butils import removeAll
    from batoms.bio import read
    removeAll()
    tio2 = read('datas/tio2_10.xyz', index = ':')
    tio2.set_frames()   
    tio2.repeat([2, 2, 2])
    tio2.model_type = 1
    assert len(tio2) == 48


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


def test_get_geometry():
    from ase.build import molecule
    from batoms.butils import removeAll
    removeAll()
    atoms = molecule('H2O')
    h2o = Batoms(atoms = atoms, label = 'h2o')
    angle = h2o.get_angle('H', 0, 'O', 0, 'H', 1)
    d = h2o.get_distances('H', 0, 'H', 1)
    com = h2o.get_center_of_mass()
    cog = h2o.get_center_of_geometry()

def test_make_real():
    from ase.build import molecule
    from batoms.butils import removeAll
    removeAll()
    atoms = molecule('H2O')
    h2o = Batoms(atoms = atoms, label = 'h2o')
    h2o.make_real()

if __name__ == '__main__':
    test_from_batom()
    test_batoms()
    test_from_coll()
    test_ase_species()
    test_get_and_set_positions()
    test_cavity()
    test_get_geometry()
    test_repeat()
    test_repeat_animation()
    print('\n Batoms: All pass! \n')