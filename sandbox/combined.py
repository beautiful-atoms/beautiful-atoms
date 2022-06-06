from batoms.batoms import Batoms
from ase.build import molecule, fcc111
au111 = fcc111("Au", (4, 4, 4), vacuum=0)
au111 = Batoms("au111", from_ase=au111)
mol = Batoms("mol", from_ase=molecule("CO"))
mol.translate([5, 5, 10])
combined = au111 + mol
combined.cell[2, 2] += 10
au111.selects["mol"].model_style = 1
au111.selects["mol"].set_scale(0.75)
au111.get_image(viewport=[1, 0, 0], output="test.png", engine="eevee")