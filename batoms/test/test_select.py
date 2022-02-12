import pytest
from batoms.butils import removeAll
from batoms.batoms import Batoms
from batoms.bio import read
import numpy as np
from time import time

def test_draw():
import numpy as np
from batoms.batoms import Batoms
from ase.build import molecule, fcc111
from batoms.butils import removeAll
removeAll()
au111 = fcc111('Au', (4, 4, 4), vacuum = 0)
mol = molecule('CH3CH2OH')
mol.translate([5, 5, 10])
atoms = au111 + mol
au111 = Batoms('au111', from_ase = atoms)
au111.cell[2, 2] += 10
indices = np.arange(64)
sel1 = au111.selects.add('ham', indices)
sel0 = au111.selects['sel0']
sel0.model_style = 1


def test():
from batoms.batoms import Batoms
from ase.build import molecule
from batoms.butils import removeAll
removeAll()
mol = molecule('H2O')
h2o = Batoms('h2o', aseAtoms = mol)
h2o.pbc = True
h2o.cell = [3, 3, 3]
h2o = h2o*[2, 1, 1]
indices = [0, 1, 2]
sel1 = h2o.selects.add('sel1', indices)
sel0 = h2o.selects['sel0']
sel0.scale = 0.5
sel1.scale = [2, 2, 2]
h2o.draw_space_filling()
sel1.model_style = 1
# h2o.split_select('sel1')
h2o.selects['sel1'].batom_dict
h2o.selects['sel1'].scale = 0.3

def test_pdb():
from batoms.batoms import Batoms
from batoms.pdbparser import read_pdb
from batoms.butils import removeAll
import numpy as np
removeAll()
atoms = read_pdb('test/datas/2piw.pdb')  # 1tim
batoms = Batoms('protein', aseAtoms = atoms)
batoms.ribbon.draw()
sel1 = batoms.selects.add('sel1', np.where(batoms.arrays['types'] == 'HETATM')[0])
sel1.show = 1
sel1.model_style = 1


if __name__ == '__main__':
    test_sheet()
    print('\n Bondsetting: All pass! \n')