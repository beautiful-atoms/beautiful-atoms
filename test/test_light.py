import pytest
import numpy as np
from batoms.utils.butils import removeAll
from batoms.render import Render, Lights, Light
from batoms.build import molecule

def test_lights():
    """
    """
    removeAll()
    r = Render(label = 'h2o')
    assert isinstance(r.lights, Lights)
    # properties
    r.lights.add('1')
    r.lights['1'].type = 'POINT'
    assert r.lights['1'].type == 'POINT'
    r.lights.add('2')
    r.lights.remove('2')
    assert len(r.lights) == 2

def test_light_direction():
    removeAll()
    h2o = molecule('h2o', 'H2O')
    h2o.render.lights['Default'].type = 'POINT'
    h2o.render.lights['Default'].energy = 1000
    assert h2o.render.lights['Default'].type == 'POINT'
    h2o.get_image([0, 0, 1], output='light-direction.png')    



if __name__ == '__main__':
    test_lights()
    test_light_direction()
    print('\n light: All pass! \n')