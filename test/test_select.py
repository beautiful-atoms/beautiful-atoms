import pytest
from batoms.butils import removeAll
from batoms.batoms import Batoms
from batoms.bio import read
import numpy as np
from time import time

def test_sheet():
from batoms.batoms import Batoms
from ase.build import molecule
from batoms.butils import removeAll
removeAll()
mol = molecule('H2O')
h2o = Batoms('h2o', atoms = mol)
h2o.pbc = True
h2o.cell = [3, 3, 3]
h2o = h2o*[2, 1, 1]
indices = [0, 1, 2]
data = {'name': 'sel1', 'indices': indices}
sel1 = h2o.selects.add('sel1', indices)
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
batoms = Batoms('protein', atoms = atoms)
sel1 = batoms.selects.add('hetatm', np.where(batoms.arrays['types'] == 'HETATM')[0])
sel1.model_style = 1
batoms.ribbon.draw()


if __name__ == '__main__':
    test_sheet()
    print('\n Bondsetting: All pass! \n')