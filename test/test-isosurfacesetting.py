import pytest
from batoms.butils import removeAll
from batoms.batoms import Batoms
from batoms.bio import read
import numpy as np
from time import time


def test_diff():
    from batoms.batoms import Batoms
    from batoms.bio import read
    removeAll()
    ag_pto_fe = read('/home/xing/ase/batoms/ag-pto-fe.cube')
    ag_pto = read('/home/xing/ase/batoms/ag-pto.cube')
    ag_pto.hide = True
    fe = read('/home/xing/ase/batoms/fe.cube')
    fe.hide = True
    volume = ag_pto_fe.isosurfacesetting.volume - ag_pto.isosurfacesetting.volume - fe.isosurfacesetting.volume
    # volume = ag_pto_fe.isosurfacesetting.volume*2
    # print(volume)
    ag_pto_fe.isosurfacesetting.volume = volume
    ag_pto_fe.isosurfacesetting[1].level = 0.008
    ag_pto_fe.isosurfacesetting[2] = {'level': -0.008, 'color': [0, 0, 1, 0.8]}
    ag_pto_fe.model_type = 1
    ag_pto_fe.draw_isosurface()
    ag_pto_fe.render.run([0, 0, 1], resolution_x = 3000, engine = 'eevee', output = 'top.png')
    ag_pto_fe.render.run([1, 0, 0], resolution_x = 3000, engine = 'eevee', output = 'side.png')


if __name__ == '__main__':
    test_diff()
    print('\n Bondsetting: All pass! \n')