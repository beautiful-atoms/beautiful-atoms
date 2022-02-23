import pytest
from batoms.utils.butils import removeAll
from batoms.pdbparser import read_pdb
from batoms.batoms import Batoms


def test_ribbon():
    removeAll()
    atoms = read_pdb('datas/1ema.pdb')  # 1ema, 1tim, 4hhb
    protein = Batoms('protein', from_ase=atoms)
    protein.ribbon.draw()
    sel1 = protein.selects.add('sel1', 'sheet A-160-A-170')
    sel1.show = True
    sel1.model_style = 1


if __name__ == '__main__':
    test_ribbon()
    print('\n Ribbon: All pass! \n')
