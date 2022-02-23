"""
modified from ASE

Module to read and write atoms in PDB file format.

"""

import warnings

import numpy as np

from ase.atoms import Atoms
from ase.cell import Cell
from ase.io.espresso import label_to_symbol


def read_atom_line(line):
    """
    Read atom line from pdb format
ATOM     25  CA  LYS A   5      48.356  17.146  -2.714  1.00  0.00           C
    """

    line = line.rstrip('\n')
    type_record = line[0:6].strip()
    name = line[12:16].strip()
    altloc = line[16]
    resname = line[17:21]
    chainid = line[21]        # Not used
    seq = line[22:26].split()
    if len(seq) == 0:
        resseq = 1
    else:
        resseq = int(seq[0])  # sequence identifier
    # icode = line[26]          # insertion code, not used
    # atomic coordinates
    try:
        coord = np.array([float(line[30:38]),
                          float(line[38:46]),
                          float(line[46:54])], dtype=np.float64)
    except ValueError:
        raise ValueError("Invalid or missing coordinate(s)")
    # occupancy & B factor
    try:
        occupancy = float(line[54:60])
    except ValueError:
        occupancy = None  # Rather than arbitrary zero or one
    if occupancy is not None and occupancy < 0:
        warnings.warn("Negative occupancy in one or more atoms")
    try:
        bfactor = float(line[60:66])
    except ValueError:
        # The PDB use a default of zero if the data is missing
        bfactor = 0.0
    # segid = line[72:76] # not used
    symbol = line[76:78].strip().upper()
    return symbol, name, altloc, resname, coord,\
        occupancy, bfactor, resseq, chainid, type_record


def read_line_cyrstal(line):
    cellpar = [float(line[6:15]),  # a
               float(line[15:24]),  # b
               float(line[24:33]),  # c
               float(line[33:40]),  # alpha
               float(line[40:47]),  # beta
               float(line[47:54])]  # gamma
    cell = Cell.new(cellpar)
    pbc = True
    return cell, pbc


def read_line_sheet(line):
    """
SHEET    1   A 9 PHE A   6  TRP A  12  0
"""
    # Chain identifier
    # sheetId = int(line[6:10].strip())
    # chainId = line[11:14].strip()
    startChain = line[21]
    # Residue sequence number
    startResi = int(line[22:26])
    endChain = line[32]
    endResi = int(line[33:37])
    sheet = {
        'name': '%s-%s-%s-%s' % (startChain, startResi, endChain, endResi),
        # 'sheetId': sheetId,
        # 'chainId': chainId,
        'startChain': startChain,
        'startResi': startResi,
        'endChain': endChain,
        'endResi': endResi,
    }
    return sheet


def read_line_helix(line):
    """
HELIX    3   3 GLY A  724  LEU A  728  5   5
    """
    startChain = line[19]
    startResi = int(line[21:25])
    endChain = line[31]
    endResi = int(line[33:37])
    helix = {
        'name': '%s-%s-%s-%s' % (startChain, startResi, endChain, endResi),
        # 'sheetId': sheetId,
        # 'chainId': chainId,
        'startChain': startChain,
        'startResi': startResi,
        'endChain': endChain,
        'endResi': endResi,
    }
    return helix


def read_pdb(fileobj, index=-1, read_arrays=True):
    """Read PDB files."""
    if isinstance(fileobj, str):
        fileobj = open(fileobj)
    images = []
    orig = np.identity(3)
    trans = np.zeros(3)
    occ = []
    bfactor = []
    types = []
    residuenames = []
    residuenumbers = []
    atomtypes = []
    chainids = []
    sheet = []
    helix = []

    symbols = []
    positions = []
    cell = None
    pbc = None

    def build_atoms():
        atoms = Atoms(symbols=symbols,
                      cell=cell, pbc=pbc,
                      positions=positions)

        if not read_arrays:
            return atoms

        info = {'occupancy': occ,
                'bfactor': bfactor,
                'residuenames': residuenames,
                'atomtypes': atomtypes,
                'residuenumbers': residuenumbers,
                'chainids': chainids,
                'types': types,
                }
        for name, array in info.items():
            if len(array) == 0:
                pass
            elif len(array) != len(atoms):
                warnings.warn('Length of {} array, {}, '
                              'different from number of atoms {}'.
                              format(name, len(array), len(atoms)))
            else:
                atoms.set_array(name, np.array(array))
        atoms.info['sheet'] = sheet
        atoms.info['helix'] = helix
        return atoms

    for line in fileobj.readlines():
        if line.startswith('CRYST1'):
            cell, pbc = read_line_cyrstal(line)
        for c in range(3):
            if line.startswith('ORIGX' + '123'[c]):
                orig[c] = [float(line[10:20]),
                           float(line[20:30]),
                           float(line[30:40])]
                trans[c] = float(line[45:55])
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # line_info = symbol, name, altloc, resname, coord, occupancy,
            #             bfactor, resseq
            line_info = read_atom_line(line)

            try:
                symbol = label_to_symbol(line_info[0])
            except (KeyError, IndexError):
                symbol = label_to_symbol(line_info[1])

            position = np.dot(orig, line_info[4]) + trans
            atomtypes.append(line_info[1])
            residuenames.append(line_info[3])
            if line_info[5] is not None:
                occ.append(line_info[5])
            bfactor.append(line_info[6])
            residuenumbers.append(line_info[7])
            chainids.append(line_info[8])
            types.append(line_info[9])

            symbols.append(symbol)
            positions.append(position)

        if line.startswith("SHEET"):
            sheet.append(read_line_sheet(line))
        if line.startswith("HELIX"):
            helix.append(read_line_helix(line))
        if line.startswith("CONECT"):
            pass
        if line.startswith("REMARK"):
            pass
        if line.startswith("HEADER"):
            pass
        if line.startswith("TITLE"):
            pass
        if line.startswith("COMPND"):
            pass
        if line.startswith('END'):
            # End of configuration reached
            # According to the latest PDB file format (v3.30),
            # this line should start with 'ENDMDL' (not 'END'),
            # but in this way PDB trajectories from e.g. CP2K
            # are supported (also VMD supports this format).
            atoms = build_atoms()
            images.append(atoms)
            occ = []
            bfactor = []
            residuenames = []
            atomtypes = []
            symbols = []
            positions = []
            cell = None
            pbc = None

    if len(images) == 0:
        atoms = build_atoms()
        images.append(atoms)
    return images[index]


if __name__ == "__main__":
    images = read_pdb('test/datas/1tim.pdb')
    # images = read_pdb('test/datas/2piw.pdb')
    print(images.arrays.keys())
    print(images.info['helix'])
