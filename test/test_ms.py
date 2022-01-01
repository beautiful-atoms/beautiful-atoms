from numpy.core.numeric import indices
import pytest
from time import time

def test_voronoi_bulk():
    from ase.build import bulk
    from batoms.batoms import Batoms
    from batoms.butils import removeAll
    removeAll()
    au = bulk('Au', cubic = True)
    au = Batoms('au', atoms = au)
    au = au*[4, 4, 4]
    au.scale = 0.1
    au.draw_voronoi()

def test_voronoi_molecule():
    from ase.build import molecule
    from batoms.batoms import Batoms
    from batoms.butils import removeAll
    removeAll()
    c2h6 = molecule('C2H6')
    c2h6.center(3.0)
    c2h6 = Batoms('c2h6', atoms = c2h6)
    c2h6.scale = 0.1
    c2h6.draw_voronoi()


def test_SAS():
    """
    isosurface sasurface 1.4

    isosurface area set 0
    36.107
    """
from ase.build import molecule
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
h2o = molecule('H2O')
h2o = Batoms('h2o', atoms = h2o)
h2o.mssetting.draw_SAS(probe = 1.4)
area = h2o.mssetting.get_sasa(partial = True)[0]

def test_SES():
    """
    isosurface sasurface 1.4

    isosurface area set 0
    36.107
    """
from ase.build import molecule
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
h2o = molecule('H2O')
h2o = Batoms('h2o', atoms = h2o)
h2o.mssetting.draw_SES(resolution=0.4, probe=1.4)
area = h2o.mssetting.get_sasa(partial = True)[0]

def perf():
tstart = time()
indices = range(1e6)
distance = np.zeros(1e6)
for i in indices:
    distance = 0
print('time: %s'(time() - tstart))


def test_SAS_area_kras():
    """
    isosurface resolution 5 sasurface 1.4

    isosurface area set 0

    Jmol    8150.5
    Batoms: 7990.7 (resolution 0.4)
    Batoms: 8132.7 (resolution 0.2)
    Batoms: 8202.8 (resolution 0.1)
    Batoms: 8232.8 (resolution 0.06)
    pymol: 8298.875
    """
from ase.io import read
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
kras = read('test/datas/kras.pdb')
kras = Batoms('kras', atoms = kras)
kras.mssetting.draw_SAS(parallel=6)
area = kras.mssetting.get_sasa()


def test_SAS():
    """
    isosurface sasurface 1.4

    isosurface area set 0
    36.107
    """
from ase.build import molecule
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
h2o = molecule('H2O')
h2o = Batoms('h2o', atoms = h2o)
h2o.draw_SAS(probe = 0.01, resolution = 0.4, threshold = 1e-6, area = True)
area = h2o.mssetting.get_sasa(partial = True)[0]

def test_SAS_area_atp():
    """

    isosurface resolution 3 sasurface 1.4

    isosurface area set 0

    Jmol    650.3
    Batoms:  655.6
    freesasa: 657.80
    """
from ase.io import read
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
atp = read('test/datas/ATP.pdb')
atp = Batoms('atp', atoms = atp)
atp.draw_SAS(probe = 1.4, resolution = 0.4, threshold = 1e-6)
area = atp.mssetting.get_sasa()
parea = atp.mssetting.get_psasa()
    assert abs(area - 657.8) < 10
    atp.get_image([0, 0, 1], padding = 3, engine = 'eevee', output = 'sas-atp.png')

def test_SAS_select():
    """

    isosurface resolution 3 sasurface 1.4

    isosurface area set 0

    Jmol    650.3
    Batoms:  655.6
    freesasa: 657.80
    """
from ase.io import read
from batoms.batoms import Batoms
from batoms.butils import removeAll
import numpy as np
removeAll()
kras = read('test/datas/kras.pdb')
kras = Batoms('kras', atoms = kras)
indices = np.where(kras.positions[:, 2] < 0)[0]
kras.draw_SAS(probe = 1.4, indices = indices, resolution = 0.4, area = True)

def test_SAS_area_kras():
    """
    isosurface resolution 5 sasurface 1.4

    isosurface area set 0

    Jmol    8150.5
    Batoms: 7990.7 (resolution 0.4)
    Batoms: 8132.7 (resolution 0.2)
    Batoms: 8202.8 (resolution 0.1)
    Batoms: 8232.8 (resolution 0.06)
    pymol: 8298.875
    """
from ase.io import read
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
kras = read('test/datas/kras.pdb')
kras = Batoms('kras', atoms = kras)
kras.draw_SAS(probe = 1.4, resolution = 0.4, area = True)
area = kras.mssetting.get_sasa()
area = kras.mssetting.get_psasa()
assert abs(area - 657.8) < 10

def test_SAS_animation():
    """

    """
from ase.io import read
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
mol = read('test/datas/deca_ala_md-pos-1.xyz', index = ':')
mol = Batoms('mol', atoms = mol, movie=True)
mol.draw_SAS(probe = 1.4, resolution = 0.4)
    areas = mol.mssetting.get_sasa()
    assert len(areas) == 10

def test_SAS_performance():
    from time import time
    from ase.io import read
    from batoms.batoms import Batoms
    from batoms.butils import removeAll
    removeAll()
    atp = read('datas/ATP.pdb')
    atp.center(5.0)
    atp = atp*[3, 3, 3]
    print(len(atp))
    atp = Batoms('atp', atoms = atp)
    tstart = time()
    atp.draw_SAS(probe = 1.4, resolution = 0.4, threshold = 1e-6, area = True)
    t = time() - tstart
    assert t < 4
    print('time: %s'%t)

def test_SES_2():
    """
    isosurface sasurface 1.4

    isosurface area set 0
    """
from ase.atoms import Atoms
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
o2 = Atoms('O2', [[-0.5, 0, 0], [0.5, 0, 0]])
o2 = Batoms('o2', atoms = o2)
o2.draw_SES(probe = 1.4, resolution = 0.4, subdivide=0)

def test_SES_3():

from ase.atoms import Atoms
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
o3 = Atoms('O3', [[-2.5, 0, 0], [2.5, 0, 0], [0, 2.5, 0]])
o3 = Batoms('o3', atoms = o3)
o3.draw_SES(probe = 1.4, resolution = 0.4, subdivide=1)


from ase.atoms import Atoms
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
ho = Atoms('HO', [[0.5, 0, 0], [-0.5, 0, 0]])
ho = Batoms('ho', atoms = ho)
# ho.draw_SAS(probe = 1.4, resolution = 0.2, threshold = 1e-6)
ho.draw_SES(probe = 1.4, resolution = 0.4, subdivide=1)


def test_SES():
    """
    isosurface sasurface 1.4

    isosurface area set 0
    """
from ase.build import molecule
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
h2o = molecule('CH3CH2OH')
print(h2o.positions)
h2o = Batoms('h2o', atoms = h2o)
# h2o.draw_SAS(probe = 1.4, resolution = 0.2, threshold = 1e-6)
h2o.draw_SES(resolution = 0.2, subdivide=0, steps = 10, probe = 1.4)

def test_SES_area_atp():
    """
    isosurface resolution 4 molecular 1.4

    isosurface area set 0
    Jmol    368.1
    Batoms: 367.3
    """
from ase.io import read
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
atp = read('test/datas/ATP.pdb')
atp = Batoms('atp', atoms = atp)
atp.draw_SES(resolution = 0.4, steps = 10, subdivide = 0)
area = atp.mssetting.get_sesa()
    assert abs(area - 368.1) < 10

def test_SES_area_kras():
    """
    isosurface resolution 3 molecular 1.4

    isosurface area set 0
            area                    volume
    Jmol    7054.5  (esolution 1)
    Jmol    17479.0  (esolution 3)
    Jmol    18518.3  (esolution 5)
    Jmol    19180.0  (esolution 8)
    Batoms: 6811   (resolution 0.4) 24134
    Batoms: 6836   (resolution 0.2) 24065
    Batoms: 7249   (resolution 0.1)  23789
    QuickSES: 6783 (resolution 0.4)
    QuickSES: 7069 (resolution 0.2)
    QuickSES: 7298 (resolution 0.1)  23211
    EDTSurf 7338 (Use slightly different vdw radii)
    """
from ase.io import read
from batoms.batoms import Batoms
from batoms.butils import removeAll
removeAll()
kras = read('test/datas/kras.pdb')
kras = Batoms('kras', atoms = kras)
kras.mssetting.draw_SES(parallel=2)
area = kras.mssetting.get_sesa()
    assert abs(area - 20016.1) < 200

def test_SES_performance():
    from time import time
    from ase.io import read
    from batoms.batoms import Batoms
    from batoms.butils import removeAll
    removeAll()
    atp = read('datas/ATP.pdb')
    atp.center(4)
    atp = atp*[3, 3, 3]
    atp = Batoms('atp', atoms = atp)
    tstart = time()
    atp.draw_SES(probe = 1.2, resolution = 0.4, threshold = 1e-4)
    t = time() - tstart
    assert t < 5
    print('time: %s'%t)

if __name__ == '__main__':
    test_voronoi_bulk()
    test_voronoi_molecule()
    test_SAS()
    test_SAS_area_atp()
    test_SAS_animation()
    test_SAS_performance()
    test_SES_area_atp()
    test_SES_performance()
    print('\n Bondsetting: All pass! \n')