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
    nframe = len(c2h6so_10.frames)
    assert nframe == 10
    # add more image
    for i in range(10, 20):
        temp = atoms.copy()
        temp.rotate(18*i, 'z')
        images.append(temp)
    c2h6so_20 = Batoms(label = 'c2h6so_20', atoms = images)
    c2h6so_20.set_frames(images)
    nframe = len(c2h6so_20.frames)
    assert nframe == 20
    write('c2h6so_10.xyz', c2h6so_10.frames)
    write('c2h6so_20.xyz', c2h6so_20.frames)
    # c2h6so.render(animation = True)

if __name__ == '__main__':
    test_batoms_animation()
    print('\n Animation: All pass! \n')