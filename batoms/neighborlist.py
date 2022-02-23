"""
"""
import numpy as np
from time import time
from ase.geometry import wrap_positions


def neighbor_kdtree(quantities, species, positions, cell, pbc,
                    cutoffs, self_interaction=False):
    """
    todo: support atoms outside cell.
    """
    tstart = time()
    maxCutoff = np.array([x for x in cutoffs.values()]).max()
    boundary = np.ones((3, 2))*maxCutoff
    boundary_data = build_boundary(species, positions, 
            cell, pbc, boundary, include_self = True)
    positions_b = boundary_data['positions']
    indices_b = boundary_data['indices']
    species_b = boundary_data['species']
    offsets_b = boundary_data['offsets']
    # print('build boundary: {:1.2f}'.format(time() - tstart1))
    #
    i = []
    j = []
    j_b = []
    # positions_b = np.dot(scaled_positions_b, cell)
    for pair, cutoff in cutoffs.items():
        indices1 = np.where(species == pair[0])[0]
        indices2 = np.where(species_b == pair[1])[0]
        p1 = positions[indices1]
        p2 = positions_b[indices2]
        if len(p1) == 0 or len(p2) == 0: continue
        indices_min = None
        # qurey the less one is faster, flip is needed
        flip = False
        if len(p1) < len(p2):
            # find max
            indices_max = primitive_neighbor_kdtree(p1, 
                p2, cutoff = cutoff[1], parallel = 1)
            # find min
            if cutoff[0] > 1e-6:
                indices_min = primitive_neighbor_kdtree(p1, 
                    p2, cutoff = cutoff[0], parallel = 1)
        else:
            flip = True
            indices_max = primitive_neighbor_kdtree(p2, 
                p1, cutoff = cutoff[1], parallel = 1)
            if cutoff[0] > 1e-6:
                indices_min = primitive_neighbor_kdtree(p2, 
                    p1, cutoff = cutoff[0], parallel = 1)
        
        n = len(p2) if flip else len(p1)
        i2 = []
        i1 = []
        if indices_min is not None:
            for k in range(n):
                indices_mid = set(indices_max[k]) - set(indices_min[k])
                m = len(indices_mid)
                i1.extend([k]*m)
                i2.extend(indices_mid)
        else:
            for k in range(n):
                # offsets1 = offsets[indices[k]]
                m = len(indices_max[k])
                i1.extend([k]*m)
                i2.extend(indices_max[k])
        # map indices to original atoms, and the boudary atoms
        if flip:
            i1_o = indices1[i2]
            i2_b = indices2[i1]
        else:
            i1_o = indices1[i1]
            i2_b = indices2[i2]
        i2_o = indices_b[i2_b]
        # remove bothways for same species, e.g. ('C', 'C')
        if pair[0] == pair[1]:
            mask = np.where((i1_o > i2_o) & ((offsets_b[i2_b] == 0).all(axis = 1)), False, True)
            i1_o = i1_o[mask]
            i2_o = i2_o[mask]
            i2_b = i2_b[mask]
        i.extend(i1_o)
        j.extend(i2_o)
        j_b.extend(i2_b)
    # Compute distance vectors.
    # print(indices3[i1][:, 0])
    tstart1 = time()
    distance_vector = positions[i] - positions_b[j_b]
    # distance_vectors.append(distance_vector)
    distances = np.sqrt(np.sum(distance_vector*distance_vector, axis = 1))
    offsets = offsets_b[j_b]
    # Remove all self-interaction.
    i = np.array(i)
    j = np.array(j)
    if not self_interaction:
        mask = np.where((i == j) & (np.prod(offsets == [0, 0, 0], axis = -1)), False, True)
        # print(mask)
        i = i[mask]
        j = j[mask]
        distances = distances[mask]
        offsets = offsets[mask]
    # print('Build distances: {:1.2f}'.format(time() - tstart))
    #=====================================
    retvals = []
    for q in quantities:
        if q == 'i':
            retvals += [i]
        elif q == 'j':
            retvals += [j]
        elif q == 'D':
            retvals += [distance_vector]
        elif q == 'd':
            retvals += [distances]
        elif q == 'S':
            retvals += [offsets]
        else:
            raise ValueError('Unsupported quantity specified.')
    print('Build bondlist: {:1.2f}'.format(time() - tstart))
    if len(retvals) == 1:
        return retvals[0]
    else:
        return tuple(retvals)

def primitive_neighbor_kdtree(positions1, positions2, 
                cutoff = 1.0, parallel = 2):
    """
    """
    from scipy.spatial import KDTree
    #
    tstart = time()
    tree = KDTree(positions2)
    indices = tree.query_ball_point(positions1, r = cutoff, workers=parallel)
    # print('KDTree: {:1.2f}'.format(time() - tstart))
    return indices



def cellPlanes(cell, origin = np.array([0, 0, 0])):
    """
    build six planes for the six faces of a cell
    """
    c1, c2, c3 = cell
    normals = np.zeros((6, 3))
    points = np.zeros((6, 3))
    normals[0] = np.cross(c2, c3)
    normals[1] = -normals[0]
    normals[2] = np.cross(c3, c1)
    normals[3] = -normals[2]
    normals[4] = np.cross(c1, c2)
    normals[5] = -normals[4]
    norm = np.linalg.norm(normals, axis=1) + 1e-8
    normals = normals/norm[:, None]
    points[0] = origin
    points[1] = origin + c1
    points[2] = origin
    points[3] = origin + c2
    points[4] = origin
    points[5] = origin + c3
    planes = []
    for i in range(3):
        plane = [{'normals': normals[2*i], 'points': points[2*i]},
                 {'normals': normals[2*i + 1], 'points': points[2*i + 1]}]
        planes.append(plane)
    return planes

def pointCellDistance(points, cell, origin = np.zeros(3)):
    """
    get distance between point and the six faces of a cell.
    """
    planes = cellPlanes(cell, origin = origin)
    distances = []
    npoint = len(points)
    distances = np.zeros((3, 2, npoint))
    for i in range(3):
        for j in range(2):
            vec = points - planes[i][j]['points']
            distances[i, j, :] = np.dot(vec, planes[i][j]['normals'])
    return distances

def build_boundary(species, positions, cell, pbc, boundary, include_self = False):
    """
    find region outside cell for periodic caculation, supercell 
    
    todo: boundary > 1
    """
    from functools import reduce
    tstart = time()
    wraped_positions = wrap_positions(positions, cell, pbc=pbc)
    distances = pointCellDistance(wraped_positions, cell)
    ib = np.array([[0, 1], [0, 1], [0, 1]])
    natom = len(positions)
    ind = np.arange(natom)
    # find atoms close to cell face with distance of r
    indices = [[], [], []]
    nb = 0
    for i in range(3):
        if not pbc[i]: continue
        ind0 =  np.where((distances[i, 0, :] < boundary[i, 1]))[0]
        ind1 =  np.where((distances[i, 1, :] < boundary[i, 0]))[0]
        indices[i] = [ind0, ind1]
        nb += len(ind0) + len(ind1)
    # print(indices)
    # init
    # repeat the scaled_positions so that it completely covers the boundary
    positions_b = np.zeros((27*nb, 3))
    offsets_b = np.zeros((27*nb, 3), dtype=int)
    indices_b = np.zeros(27*nb, dtype=int)
    species_b = np.zeros(27*nb, dtype = 'U20')
    # build face
    offset = [1, -1]
    nb1 = 0
    nb2 = 0
    for c in range(3):
        if not pbc[i]: continue
        for i in range(2):
            nb2 = nb1 + len(indices[c][i])
            indices_b[nb1:nb2] = indices[c][i]
            positions_b[nb1:nb2] = positions[indices[c][i]]
            offsets_b[nb1:nb2][:, c] = offset[i]
            nb1 = nb2
    # build edge
    for c in [[0, 1], [0, 2], [1, 2]]:
        if not pbc[c[0]] and not pbc[c[1]]: continue
        for i in range(2):
            for j in range(2):
                indices2 = np.intersect1d(indices[c[0]][i], indices[c[1]][j])
                nb2 = nb1 + len(indices2)
                indices_b[nb1:nb2] = indices2
                positions_b[nb1:nb2] = positions[indices2]
                offsets_b[nb1:nb2][:, c[0]] = offset[i]
                offsets_b[nb1:nb2][:, c[1]] = offset[j]
                nb1 = nb2
    # build corner
    if pbc[0] and pbc[1] and pbc[2]:
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    indices3 = reduce(np.intersect1d, 
                            (indices[0][i], indices[1][j], indices[2][k]))
                    nb2 = nb1 + len(indices3)
                    indices_b[nb1:nb2] = indices3
                    positions_b[nb1:nb2] = positions[indices3]
                    offsets_b[nb1:nb2][:, 0] = offset[i]
                    offsets_b[nb1:nb2][:, 1] = offset[j]
                    offsets_b[nb1:nb2][:, 2] = offset[k]
                    nb1 = nb2
    positions_b = positions_b[0:nb2]
    indices_b = indices_b[0:nb2]
    species_b = species[indices_b]
    offsets_b = offsets_b[0:nb2]
    positions_b = positions_b + np.dot(offsets_b, cell)
    if include_self:
            positions_b = np.append(positions, positions_b, axis = 0)
            indices_b = np.append(np.arange(natom), indices_b)
            species_b = np.append(species, species_b)
            offsets_b = np.append(np.zeros((natom, 3)), offsets_b, axis = 0)

    # print('build boundary: {:1.2f}'.format(time() - tstart))
    boundary_data = {
        'positions': positions_b,
        'indices': indices_b,
        'species': species_b,
        'offsets': offsets_b,
    }
    return boundary_data

def test_kdtree(atoms, cutoffs1, cutoffs2):
    
    tstart = time()
    i, j, d, S = neighbor_list('ijdS', atoms, cutoffs1)
    print('time %s'%(time() - tstart))
    if len(i) < 20:
        print(i, j, d)
    tstart = time()
    species = atoms.numbers
    i, j, d, offsets = neighbor_kdtree('ijdS', species, 
                atoms.positions, atoms.get_cell(complete=True), atoms.pbc,
            cutoffs2)
    print('time %s'%(time() - tstart))
    if len(i) < 20:
        print(i, j, d)

def test_pointCellDistance():
    from ase.io import read, write
    from ase.build import molecule, bulk
    from ase.visualize import view
    atoms = molecule('H2O')
    # atoms.center(3)
    # atoms = bulk('Pt', cubic=True)
    atoms.pbc = True
    positions = atoms.positions
    distances = pointCellDistance(positions, atoms.cell)
    print(distances)

def test_buildBoundary(atoms, cutoffs):
    maxCutoff = np.array([x for x in cutoffs.values()]).max()
    positions = atoms.positions
    species = np.array(atoms.get_chemical_symbols())
    positions_b, indices_b, species_b, offsets_b =  \
            build_boundary(species, positions, atoms.cell, atoms.pbc, maxCutoff)
    face = Atoms(species_b, positions_b)
    # face = Atoms(['O']*len(species_b), positions_b + np.dot(offsets_b, atoms.cell))
    # print(face)
    # view(atoms + face)

if __name__ == "__main__":
    from ase.io import read, write
    from ase.build import molecule, bulk
    from ase.visualize import view
    from ase import Atom, Atoms
    from ase.neighborlist import neighbor_list
    atoms = molecule('H2O')
    # atoms.pbc = True
    # atoms.center(0.5)
    # atoms = bulk('Au', cubic = True)
    # atoms.cell = [4, 7, 5]
    # atoms = read('test/datas/mof-5.cif')
    # atoms = atoms*[4, 4, 4]
    # atoms = read('test/datas/1tim.pdb')
    print(len(atoms))
    # view(atoms)
    cutoffs1 = {
        ('O', 'H'): 1.261, 
        ('C', 'H'): 1.261, 
        ('C', 'C'): 1.261, 
            ('Au', 'Au'):3.0,
                }
    cutoffs2 = {
        (8, 1)  : [0.000, 1.261],
        (6, 1)  : [0.000, 1.261],
        (6, 6)  : [0.000, 1.261],
        (79, 79): [0, 3.0],
        # (1, 8)  : [1.200, 2.100],
        }
    # test_buildBoundary(atoms, cutoffs2)
    test_kdtree(atoms, cutoffs1, cutoffs2)