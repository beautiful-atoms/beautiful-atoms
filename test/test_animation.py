
def test_animation():
from ase.io import read, write
from batoms import Batoms
from batoms.butils import removeAll
removeAll()
atoms = read('test/datas/deca_ala_md-pos-1.xyz', index = ':')
batoms = Batoms('c2h6so', from_ase = atoms, movie = True)
batoms.model_style = 1

def test_animation_2():
    from ase.io import read, write
    from batoms import Batoms
    from batoms.butils import removeAll
    removeAll()
    atoms = read('datas/tio2_10.xyz', index = ':')
    tio2 = Batoms('tio2', atoms = atoms, movie = True)
    tio2.boundary = 0.01
    tio2.model_style = 1

def test_animation_metaball():
from ase.io import read
from batoms import Batoms
from batoms.butils import removeAll
removeAll()
atoms = read('test/datas/tio2_10.xyz', index = ':')
batoms = Batoms('c2h6so', atoms = atoms, movie = True, shape = 3)
    # batoms.model_style = 1


def test_batoms_animation():
    from ase.build import molecule
    from batoms.butils import removeAll
    from batoms import Batoms
    from ase.io import write
    removeAll()
    atoms = molecule('C2H6SO')
    images = []
    for i in range(10):
        temp = atoms.copy()
        temp.rotate(18*i, 'z')
        images.append(temp)
    
    c2h6so_10 = Batoms(label = 'c2h6so_10', atoms = images)
    c2h6so_10.set_frames(images)
    #
    assert c2h6so_10.nframe == 10
    # add more image
    for i in range(10, 20):
        temp = atoms.copy()
        temp.rotate(18*i, 'z')
        images.append(temp)
    c2h6so_20 = Batoms(label = 'c2h6so_20', atoms = images)
    c2h6so_20.set_frames(images)
    assert c2h6so_20.nframe == 20
    write('datas/c2h6so_10.xyz', c2h6so_10.frames)
    # c2h6so.render(animation = True)


def test_batoms_animation_2():
    from batoms.butils import removeAll
    from batoms import Batoms
    from ase.io import write, read
    removeAll()
    atoms = read('datas/tio2.cif')
    atoms = atoms*[2, 2, 2]
    images = []
    for i in range(100):
        temp = atoms.copy()
        temp[0].x += i*0.01
        images.append(temp)
    
    tio2 = Batoms(label = 'tio2', atoms = images)
    tio2.set_frames(images)
    tio2.boundary = 0.01
    tio2.model_style = 1
    write('datas/tio2_30.xyz', tio2.frames)


if __name__ == '__main__':
    test_animation()
    test_animation_2()
    test_batoms_animation()
    test_batoms_animation_2
    print('\n Animation: All pass! \n')