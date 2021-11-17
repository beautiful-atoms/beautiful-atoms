import pytest
import numpy as np

def test_light():
    """
    """
    from batoms.butils import removeAll
    from batoms.light import Lights, Light
    removeAll()
    light = Light(label = 'h2o', name = '1')
    assert isinstance(light, Light)
    # properties
    light.type = 'POINT'
    assert light.type == 'POINT'
    light.energy = 10
    assert light.energy == 10
def test_lights():
    """
    """
    from batoms.butils import removeAll
    from batoms.light import Lights, Light
    removeAll()
    lights = Lights(label = 'h2o')
    assert isinstance(lights, Lights)
    # properties
    lights.add('1')
    lights['1'].type = 'POINT'
    assert lights['1'].type == 'POINT'
    lights.add('2')
    lights.remove('2')
    assert len(lights) == 1

def test_render():
    from batoms.build import molecule
    from batoms.butils import removeAll
    removeAll()
    h2o = molecule('h2o', 'H2O')
    h2o.render.lights.add('right', direction = [1, 0, 0])
    assert len(h2o.render.lights) == 2
    h2o.render.lights['right'].type = 'POINT'
    assert h2o.render.lights['right'].type == 'POINT'
    h2o.render.lights['Default'].direction = [1, 0, 0]
    h2o.render.run([0, 0, 1], resolution_x = 200, output='light-direction.png')    



if __name__ == '__main__':
    test_light()
    test_lights()
    test_render()
    print('\n light: All pass! \n')