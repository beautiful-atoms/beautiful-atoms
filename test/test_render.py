import pytest
import numpy as np



def test_render():
    """
    """
    from batoms.render import Render
    from batoms.butils import removeAll
    removeAll()
    render = Render('test')
    assert isinstance(render, Render)
    render.studiolight = 'paint.sl'
    render.viewport = [1, 0, 0]
    render.resolution = [1000, 1000]
    # from collection
    render = Render('test')
    assert isinstance(render, Render)
    assert render.studiolight == 'paint.sl'
    assert render.resolution == [1000, 1000]


def test_render_init():
    from batoms.butils import removeAll
    from batoms import Batoms
    removeAll()
    from ase.build import molecule
    h2o = molecule('H2O')
    h2o = Batoms('h2o', atoms = h2o)
    h2o.render.init()
    h2o.render.viewport = [1, 0, 0]
    nh3 = molecule('NH3')
    nh3 = Batoms('nh3', atoms = nh3)
    nh3.translate([3, 3, 0])
    nh3.render.init()
    nh3.render.viewport = [0, 1, 0]
    nh3.render.lights.add('1', type = 'SUN', direction = [1, 0, 0])
    nh3.get_image(output='render-init.png')

def test_render_setter():
    from batoms.butils import removeAll
    from batoms import Batoms
    from batoms.render import Render
    removeAll()
    from ase.build import molecule
    h2o = molecule('H2O')
    h2o = Batoms('h2o', atoms = h2o)
    h2o.render.init()
    h2o.render.viewport = [1, 0, 0]
    nh3 = molecule('NH3')
    nh3 = Batoms('nh3', atoms = nh3)
    nh3.translate([3, 3, 0])
    nh3.render = Render('nh3')
    nh3.render.viewport = [0, 1, 0]
    nh3.render.lights.add('1', type = 'SUN', direction = [1, 0, 0])
    nh3.get_image(output='render-setter.png')

def test_render_au111():
    from ase.build import fcc111
    from batoms.butils import removeAll
    from batoms import Batoms
    removeAll()
    atoms = fcc111('Au', size = (4, 4, 4), vacuum=0)
    au111 = Batoms(label = 'au111', atoms = atoms)
    au111.cell[2, 2] += 5
    au111.draw_cell()
    au111.get_image()
    # au111.render.run([1, 1, 1], engine = 'eevee', resolution_x = 200, output='au111-eevee')
    # au111.render.run([1, 1, 1], engine = 'cycles', resolution_x = 200, output='au111-cycles')

if __name__ == '__main__':
    test_render()
    test_render_setter()
    test_render_au111()
    print('\n render: All pass! \n')
