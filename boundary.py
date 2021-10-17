from ase import Atoms
import numpy as np
from time import time


def search_boundary(atoms, boundary = [[0, 1], [0, 1], [0, 1]], skin = 3):
    """
    search atoms in the boundary

    Parameters:

    atoms: ASE Atoms object

    boundary: list

    skin: float
        Could be the maximum cutoff.

    Return:

    atoms_boundary: ASE Atoms object
        atoms inside the boudary, not include the core
    offsets_skin: numpy array
        atoms close to boundary with a skin distance. 
        Used to search bond outside the boundary.
        index and offset

    """
    from ase.cell import Cell
    from math import floor, ceil
    tstart = time()
    if isinstance(boundary, float):
        boundary = [[-boundary, 1 + boundary], [-boundary, 1+boundary], [-boundary, 1+boundary]]
    boundary = np.array(boundary)
    boundary_skin = boundary.copy()
    # skin to scaled distance in cell
    par = atoms.cell.cellpar()
    skin = np.array([skin/par[0], skin/par[1], skin/par[2]])
    # skin region
    boundary_skin[:, 0] += skin
    boundary_skin[:, 1] -= skin
    # find supercell
    f = np.floor(boundary)
    c = np.ceil(boundary)
    ib = np.array([f[:, 0], c[:, 1]]).astype(int)
    M = np.product(ib[1] - ib[0] + 1)
    positions = atoms.get_scaled_positions()
    n = len(positions)
    npositions = np.tile(positions, (M - 1,) + (1,) * (len(positions.shape) - 1))
    i0 = 0
    # index
    offsets = np.zeros((M*n, 4), dtype=int)
    ind0 = np.arange(n).reshape(-1, 1)
    symbols0 = atoms.get_chemical_symbols()
    species0 = atoms.info['species']
    symbols = []
    species = []
    # repeat the positions so that
    # it completely covers the boundary
    for m0 in range(ib[0, 0], ib[1, 0] + 1):
        for m1 in range(ib[0, 1], ib[1, 1] + 1):
            for m2 in range(ib[0, 2], ib[1, 2] + 1):
                if m0 == 0 and m1 == 0 and m2 == 0: continue
                i1 = i0 + n
                npositions[i0:i1] += (m0, m1, m2)
                offsets[i0:i1] = np.append(ind0, np.array([[m0, m1, m2]]*n), axis = 1)
                symbols.extend(symbols0)
                species.extend(species0)
                i0 = i1
    # boundary condition
    ind1 =  np.where((npositions[:, 0] > boundary[0][0]) & (npositions[:, 0] < boundary[0][1]) \
                & (npositions[:, 1] > boundary[1][0]) & (npositions[:, 1] < boundary[1][1]) \
                & (npositions[:, 2] > boundary[2][0]) & (npositions[:, 2] < boundary[2][1]))
    # build atoms inside the boundary
    ind1 = list(set(ind1[0]))
    npositions1 = npositions[ind1]
    npositions1 = np.dot(npositions1, atoms.cell)
    symbols = np.array(symbols)[ind1]
    species = np.array(species)[ind1]
    atoms_boundary = Atoms(symbols, npositions1)
    atoms_boundary.info['species'] = species
    # build atoms inside the skin region
    # could be the atoms from core, thus add core
    npositions = np.append(npositions, positions, axis = 0)
    ind2 = np.append(ind1, i0 + ind0).astype(int)
    offsets[i0:i0+n] = np.append(ind0, np.array([[0, 0, 0]]*n), axis = 1)
    # atoms not belong to the skin region
    ind3 =  np.where((npositions[:, 0] > boundary_skin[0][0]) & (npositions[:, 0] < boundary_skin[0][1]) \
                & (npositions[:, 1] > boundary_skin[1][0]) & (npositions[:, 1] < boundary_skin[1][1])  #\
                & (npositions[:, 2] > boundary_skin[2][0]) & (npositions[:, 2] < boundary_skin[2][1]))
    ind3 = list(set(ind3[0]))
    ind3 = list(set(ind2) - set(ind3))
    offsets_skin = offsets[ind3]
    # print('search boundary: {0:10.2f} s'.format(time() - tstart))
    return atoms_boundary, offsets_skin

def search_bond(positions0, offsets_skin, bondlists, boundary, recursive = False, previous = None):
    """
    Search bonded atoms of sp1 or sp2 recursively.

    Parameters:

    positions0: array
        positions of original atoms, the core atoms
    offsets_skin: array
        index and offset of atoms of sp1 or sp2
    bondlists: array
        bondlist of the core atoms
    boundary: array
    recursive: bool
        recursive or not, for search mode 1 or 2
    previous: array
        offsets_skin of previous search. 
        The index should be remove from the result of this search.
    
    Return:

    offset_new: array
        index and offset of atoms which bond to sp1 or sp2
    
    """
    neighbors = np.array([]).reshape(-1, 4)
    if len(offsets_skin) == 0: return neighbors
    # build bonded atoms using offset_skin and bondlists
    sites = set(offsets_skin[:, 0])
    for i in sites:
        ind = np.where(bondlists[:, 0] == i)[0]
        if len(ind) == 0: continue
        neighbor = bondlists[ind][:, 1:]
        sitei = np.where(offsets_skin[:, 0] == i)[0]
        for j in sitei:
            temp = neighbor.copy()
            temp[:, 1:4] = temp[:, 1:4] + offsets_skin[j][1:4]
            neighbors = np.append(neighbors, temp, axis = 0)
    # choose index by two pricinple
    # 1, outside the boundary
    # 2, not from preious index
    # rebuild the positions, and check boundary
    neighbors = neighbors.astype(int)
    npositions = positions0[neighbors[:, 0]] + neighbors[:, 1:4]
    index =  np.where((npositions[:, 0] > boundary[0][0]) & (npositions[:, 0] < boundary[0][1]) \
                & (npositions[:, 1] > boundary[1][0]) & (npositions[:, 1] < boundary[1][1]) \
                & (npositions[:, 2] > boundary[2][0]) & (npositions[:, 2] < boundary[2][1]))
    #
    mask = np.ones(npositions.shape[0], dtype=bool)
    mask[index] = False
    offset_new = neighbors[mask]
    # offset_new = np.unique(offset_new, axis=0)
    # remove previous index, and remove duplicate
    if previous is not None:
        __, indices = np.unique(np.concatenate([previous, offset_new]), return_index=True, axis=0)
        indices = indices[indices >= len(previous)] - len(previous)
        offset_new = offset_new[indices]
    if len(offset_new) == 0: return np.array([]).reshape(-1, 4)
    # recursive search
    if recursive:
        offset_new1 = search_bond(positions0, offset_new, bondlists, boundary, recursive = True, previous = offsets_skin)
        # save
        if offset_new1 is not None:
            offset_new = np.append(offset_new, offset_new1, axis = 0)
    else:
        return offset_new
    return offset_new


class Boundary:
    """
    Boundary object.
    Example:
      cl = Boundary(atoms, d, index)
      cl.build()
    """

    def __init__(self, atoms, boundary_list, rotate_atoms = False):
        self.atoms = atoms
        self.natoms = len(atoms)
        self.cell = atoms.cell
        self.boundary_list = boundary_list
        # self.d = d
        # self.index = index
        self.rotate_atoms = rotate_atoms
        #
    def build(self, ):
        for boundary in self.boundary_list:
            print(boundary)
            self.cut(**boundary)
        return self.atoms
    def cut(self, atoms = None, d = None, index = None, direction = 1):
        """
        """
        if not atoms:
            atoms = self.atoms
        cell = atoms.cell
        normal = self.get_plane(d, index, cell)
        print(normal, d)
        # a*x + b*y + c*z - d = 0
        mask = []
        natoms = len(atoms)
        for i in range(natoms):
            v = atoms[i].position.dot(normal) - d
            # print(v)
            if v*direction > 0:
                mask.append(i)
        # view(atoms[mask])
        del atoms[mask]
        self.atoms = atoms
        if self.rotate_atoms:
            atoms = self.rotate()
    def rotate(self, atoms = None, index = None):
        """
        rotate normal of plane to z axis
        """
        import scipy
        if not atoms:
            atoms = self.atoms
        if not index:
            index = self.index
        cell = atoms.cell
        normal, d = self.get_plane(self.d, self.index, cell)
        vec = np.cross([0.0000014159, 0.000001951, 1], index)
        vec = vec/np.linalg.norm(vec)
        ang = np.arccos(normal[2])
        vec = ang*vec
        r = scipy.spatial.transform.Rotation.from_rotvec(vec)
        if scipy.version.version >= '1.4':
            mat = r.as_matrix()
        else: 
            mat = r.as_dcm()
        # print(mat)
        atoms.positions = atoms.positions.dot(mat)
        # atoms.cell = atoms.cell.dot(mat)
        self.atoms = atoms
        return atoms
    def get_plane(self, d, index = None, cell = None):
        '''
        plane equation: three point and distance from origin
        return normal and distance
        # a*x + b*y + c*z - d = 0
        '''
        index = [1.0/(index[i] + 0.000001) for i in range(3)]
        index = np.array(index)
        index = index/np.linalg.norm(index)
        points = cell*index
        # print(points)
        # x 
        v1 = points[1] - points[0]
        v2 = points[2] - points[0]
        normal = np.cross(v1, v2)
        print(v1, v2, normal)
        normal = normal/np.linalg.norm(normal)
        a, b, c = normal
        # d = np.dot(normal, points[2])
        return normal    
def build_bondlists(atoms, cutoff):
    """
    The default bonds are stored in 'default_bonds'
    Get all pairs of bonding atoms
    remove_bonds
    """
    from neighborlist import neighbor_list
    if len(cutoff) == 0: return {}
    #
    tstart = time()
    nli, nlj, nlS = neighbor_list('ijS', atoms, cutoff, self_interaction=False)
    bondlists = np.append(np.array([nli, nlj], dtype=int).T, np.array(nlS, dtype=int), axis = 1)
    # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
    return bondlists


if __name__ == "__main__":
    from ase.io import read, write
    from ase.build import molecule, bulk
    from ase.visualize import view
    from ase import Atom, Atoms
    # atoms = bulk('Pt', cubic = True)
    # atoms = atoms*[6, 6, 6]
    # cl = Boundary(atoms, boundary_list = [{'d': 10.0, 'index': [2, 2, 1]}])
    # atoms = cl.build()
    # view(atoms)
    # atoms = bulk('Pt', cubic = True)
    # atoms.write('pt.in')
    # view(atoms)
    # atoms = atoms*[3, 3, 3]
    atoms = read('docs/source/_static/datas/tio2.cif')
    atoms.info['species'] = atoms.get_chemical_symbols()
    # atoms = read('docs/source/_static/datas/mof-5.cif')
    # atoms.positions -= atoms.get_center_of_mass()
    # positions1, positions2, offset2 = search_boundary(atoms.positions, atoms.cell, boundary=[[0., 1.], [0., 1.], [0., 1.]])
    #positions1, positions2, offset2 = search_boundary(atoms.positions, atoms.cell, boundary=[[-2.6, 1.6], [-2.6, 1.6], [-2.6, 1.6]])
    cutoff = {('Ti', 'O'): [0, 3.0]}
    bondlists = build_bondlists(atoms, cutoff)
    print(bondlists)
    skins = np.array([[1, 4, 1, 1]])
    bondsetting = [('Ti', 'O')]
    boundary = np.array([[-0.6, 1.6], [-0.6, 1.6], [-0.6, 1.6]])
    # bondlists, only search bond
    # skins, only search atoms
    skins = search_bond(atoms.positions, skins, bondlists, boundary)
    print(skins)
    # skin = Atoms('Au'*len(positions2), positions = positions2)
    # view(atoms + skin)
    # bondsetting = {('Ti', 'O'): [0, 2.5, False, False], ('O', 'O'): [0, 1.5, False, False]}
    # bondpairs_boundary = search_skin(atoms, bondsetting, boundary=[[-0.6, 1.6], [-0.6, 1.6], [-0.6, 1.6]])
    # print(bondpairs_boundary)
