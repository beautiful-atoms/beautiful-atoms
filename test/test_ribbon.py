import pytest
from batoms.butils import removeAll
from batoms.batoms import Batoms
from batoms.bio import read
import numpy as np
from time import time

def test_sheet():
from batoms.batoms import Batoms
from batoms.pdbparser import read_pdb
from batoms.butils import removeAll
removeAll()
atoms = read_pdb('test/datas/2piw.pdb')  # 1tim
batoms = Batoms('protein', atoms = atoms)
batoms.ribbon.draw_sheet()
batoms.ribbon.draw_helix()
batoms.ribbon.draw_turn()



def test_helix():
from batoms.batoms import Batoms
from batoms.pdbparser import read_pdb
from batoms.butils import removeAll
removeAll()
atoms = read_pdb('test/datas/2piw.pdb')  # 1tim
batoms = Batoms('protein', atoms = atoms)
batoms.ribbon.helix.draw_helix()

if __name__ == '__main__':
    test_sheet()
    print('\n Bondsetting: All pass! \n')