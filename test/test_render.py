import pytest
import numpy as np
from ase.build import molecule, fcc111
from batoms.render import Render
from batoms import Batoms
from batoms.utils.butils import removeAll


def test_render():
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
    removeAll()
    h2o = molecule('H2O')
    h2o = Batoms('h2o', from_ase=h2o)
    h2o.render.init()
    h2o.render.viewport = [1, 0, 0]
    h2o.get_image()


def test_render_setter():
    removeAll()
    h2o = molecule('H2O')
    h2o = Batoms('h2o', from_ase=h2o)
    h2o.render.init()
    h2o.render.viewport = [1, 0, 0]
    nh3 = molecule('NH3')
    nh3 = Batoms('nh3', from_ase=nh3)
    nh3.translate([3, 3, 0])
    nh3.render = Render('nh3')
    nh3.render.viewport = [0, 1, 0]
    nh3.render.lights.add('1', type='SUN', direction=[1, 0, 0])
    nh3.get_image(output='render-setter.png')


def test_render_au111():
    removeAll()
    atoms = fcc111('Au', size=(4, 4, 4), vacuum=0)
    au111 = Batoms(label='au111', from_ase=atoms)
    au111.cell[2, 2] += 5
    au111.get_image()
    au111.get_image([1, 1, 1], engine='eevee', output='au111-eevee.png')
    au111.get_image([1, 1, 1], engine='cycles', output='au111-cycles.png')


if __name__ == '__main__':
    test_render()
    test_render_setter()
    test_render_au111()
    print('\n render: All pass! \n')
