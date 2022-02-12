import pytest
from batoms.butils import removeAll
from batoms.batoms import Batoms
from batoms.bio import read
import numpy as np
from time import time

def test_slice():
from batoms.batoms import Batoms
from batoms.bio import read
from batoms.butils import removeAll
removeAll()
h2o = read('/home/xing/ase/batoms/h2o-homo.cube')
# h2o = read('/home/xing/ase/batoms/fe.cube')
h2o.isosurfacesetting['1'] = {'level':-0.001}
h2o.isosurfacesetting['2'] = {'level':0.001, 'color': [0, 0, 0.8, 0.5]}
h2o.isosurfacesetting.draw_isosurface()
h2o.planesetting[(1, 0, 0)] = {'distance': 6, 'slicing': True}
h2o.planesetting.draw_lattice_plane()
h2o.get_image([0, 0, 1], engine = 'eevee')

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
    ag_pto_fe.model_style = 1
    ag_pto_fe.draw_isosurface()
    ag_pto_fe.render.resolution = [3000, 3000]
    ag_pto_fe.get_image([0, 0, 1], engine = 'eevee', output = 'top.png')
    ag_pto_fe.get_image([1, 0, 0], engine = 'eevee', output = 'side.png')


if __name__ == '__main__':
    test_slice()
    test_diff()
    print('\n Bondsetting: All pass! \n')