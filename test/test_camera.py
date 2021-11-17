import pytest
import numpy as np

def test_camera():
    """
    """
    from batoms.butils import removeAll
    from batoms.camera import Camera
    removeAll()
    camera = Camera(label = 'h2o', name = '1')
    assert isinstance(camera, Camera)
    # properties
    camera.type = 'PERSP'
    assert camera.type == 'PERSP'
    camera.lens = 50
    assert camera.lens == 50
    camera.target = [5, 0, 0]


def test_render():
    from batoms.build import molecule
    from batoms.butils import removeAll
    removeAll()
    h2o = molecule('h2o', 'H2O')
    h2o.render.run([1, 0, 0], output='camera.png')    
    # target
    h2o.render.camera.target = h2o['H'].positions[1]
    h2o.render.run(output = 'camera-target.png')


if __name__ == '__main__':
    test_camera()
    test_render()
    print('\n camera: All pass! \n')