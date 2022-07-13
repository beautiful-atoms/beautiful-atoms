import os
from ase import io
from ase.io.cube import read_cube_data
from batoms import Batoms
# from time import time


def read(filename, label = None, **kwargs):
    """
    wrapper function for ase.io.read
    """
    base = os.path.basename(filename)
    base = os.path.splitext(base)
    if label is None:
        label = base[0]
        label = label.replace('-', '_')
        if label[:-1].isdigit():
            label = 'b_' + label
    ext = base[1]
    if ext == '.cube':
        # tstart = time()
        volume, atoms = read_cube_data(filename, **kwargs)
        # print('Read cube: {0:1.2f}'.format(time() - tstart))
        batoms = Batoms(label, from_ase=atoms, volume={label: volume})
    elif "CHGCAR" in filename.upper():
        from pymatgen.io.vasp.outputs import VolumetricData
        poscar, data, data_aug = VolumetricData.parse_file('CHGCAR')
        # load structure and vlumetric data into Batoms
        batoms = Batoms('batoms', from_pymatgen = poscar.structure)
        batoms.volumetric_data['chgcar'] = data['total']
    else:
        atoms = io.read(filename=filename, **kwargs)
        batoms = Batoms(label=label, from_ase=atoms)
    return batoms
