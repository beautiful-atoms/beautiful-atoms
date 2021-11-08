import pytest
import numpy as np



def test_render():
    """
    """
    from batoms import Batoms
    from batoms.butils import removeAll
    from batoms.render import Render
    removeAll()
    r = Render('test')
    assert isinstance(r, Render)
    #
    from ase.build import fcc111
    removeAll()
    atoms = fcc111('Au', size = (4, 4, 4), vacuum=0)
    au111 = Batoms(label = 'au111', atoms = atoms)
    au111.cell[2, 2] += 5
    au111.draw_cell()
    au111.render.run([1, 1, 1], engine = 'workbench', resolution_x = 200, output='au111')
    au111.render.run([1, 1, 1], engine = 'eevee', resolution_x = 200, output='au111-eevee')
    au111.render.run([1, 1, 1], engine = 'cycles', resolution_x = 200, output='au111-cycles')

if __name__ == '__main__':
    test_render()
    print('\n render: All pass! \n')