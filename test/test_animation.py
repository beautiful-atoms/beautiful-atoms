from ase.io import read
from batoms import Batoms
from batoms.utils.butils import removeAll


def test_animation_molecule():
    removeAll()
    deca = read("datas/deca_ala_md.xyz", index=":")
    deca = Batoms("deca", from_ase=deca, movie=True)
    deca.model_style = 1


def test_animation_crystal():
    removeAll()
    atoms = read("datas/tio2_10.xyz", index=":")
    tio2 = Batoms("tio2", from_ase=atoms, movie=True)
    tio2.model_style = 2
    tio2.boundary = 0.01


if __name__ == "__main__":
    test_animation_molecule()
    test_animation_crystal()
    print("\n Animation: All pass! \n")
