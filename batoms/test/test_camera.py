from batoms import Batoms
import pytest
import numpy as np

def test_camera():
    """
    """
    from batoms.butils import removeAll
    from batoms.render import Camera
    removeAll()
    camera = Camera(label = 'h2o', name = '1')
    assert isinstance(camera, Camera)
    # properties
    camera.type = 'PERSP'
    assert camera.type == 'PERSP'
    camera.lens = 50
    assert camera.lens == 50
    camera.target = [5, 0, 0]


def test_camera_ortho():
    from ase.build import fcc111
    from batoms import Batoms
    from batoms.butils import removeAll
    removeAll()
    au111 = fcc111('Au', (4, 4, 4), vacuum = 0)
    au111 = Batoms('au111', atoms = au111)
    au111.cell[2, 2] += 10
    au111.draw_cell()
    au111.get_image([1, 0, 0], output='camera_ortho.png')
    # 
    au111.get_image([0, 0, 1], center = au111['Au'].positions[-1], 
                output = 'camera-center.png')
    #

def test_camera_persp():
    from ase.build import fcc111
    from batoms import Batoms
    from batoms.butils import removeAll
    removeAll()
    au111 = fcc111('Au', (4, 10, 4), vacuum = 0)
    au111 = Batoms('au111', atoms = au111)
    au111.cell[2, 2] += 10
    au111.draw_cell()
    au111.render.camera.type = 'PERSP'
    au111.get_image([1, 0, 0], output = 'camera-persp.png')
    

if __name__ == '__main__':
    test_camera()
    test_camera_ortho()
    test_camera_persp()
    print('\n camera: All pass! \n')