import pytest
from batoms.butils import removeAll
from batoms.batoms import Batoms
from batoms.bio import read
import numpy as np
from time import time

def test_select():
from batoms.batoms import Batoms
from batoms.pdbparser import read_pdb
from batoms.butils import removeAll
removeAll()
atoms = read_pdb('test/datas/1ema.pdb')  # 1ema, 1tim, 4hhb
protein = Batoms('protein', from_ase = atoms)
protein.ribbon.draw()
sel1 = protein.selects.add('sel1', 'sheet A-160-A-170')
# sel1 = batoms.selects.add('sel1', 'sheet A-148-A-152')
# sel1 = batoms.selects.add('sel1', 'helix A-94-A-112')
sel1.show = True
sel1.model_style = 1
# sel2 = protein.selects.add('sel2', 'hetatm')
# sel2.show = True


def test_sheet():
from batoms.batoms import Batoms
from batoms.pdbparser import read_pdb
from batoms.butils import removeAll
removeAll()
atoms = read_pdb('test/datas/1ema.pdb')  # 1ema, 1tim
batoms = Batoms('protein', from_ase = atoms)
batoms.ribbon.draw_sheet()
batoms.ribbon.draw_helix()
batoms.ribbon.draw_turn()
# batoms.ribbon.draw()

def test_helix():
from batoms.batoms import Batoms
from batoms.pdbparser import read_pdb
from batoms.butils import removeAll
removeAll()
atoms = read_pdb('test/datas/1jj2.pdb')  # 1tim
batoms = Batoms('protein', from_ase = atoms)
batoms.ribbon.draw()

if __name__ == '__main__':
    test_sheet()
    print('\n Bondsetting: All pass! \n')