"""
"""
import numpy as np
from time import time
from ase.geometry import wrap_positions
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)

def RemovePbc(species0, positions0, cell, pbc, cutoffs):
    """
    convert atoms with pbc to non-pbc by adding boundary
    atoms within max cutoff
    """
    maxCutoff = np.array([x for x in cutoffs.values()]).max()
    boundary = np.ones((3, 2))*maxCutoff
    array = build_boundary(species0, positions0,
                           cell, pbc, boundary, include_self=True)
    return array


def bondlist_kdtree(quantities, species0, positions0, cell, pbc,
                    setting, self_interaction=False):
    """
    return

    i: index1
    j: index2
    k: search bond style
    """
    cutoffs = {}
    for pair, b in setting.items():
        cutoffs[pair] = [b['min'], b['max']]
    natom = len(positions0)
    # orignal atoms
    array1 = {
        'positions': positions0,
        'species': species0,
        'indices': np.arange(natom),
        'offsets': np.zeros((natom, 3)),
    }
    # atoms added with boundary
    array2 = RemovePbc(species0, positions0, cell, pbc, cutoffs)
    bonddatas = primitive_neighbor_kdtree(array1, array2, cutoffs)
    # build bondlist
    i = []
    j = []
    k = []
    p = []
    t = []
    for pair, data in bonddatas.items():
        if len(data) == 0:
            continue
        i1 = []
        j1 = []
        for i2, j2 in data.items():
            n = len(j2)
            i1.extend([i2]*n)
            j1.extend(j2)
        # remove bothways for same species, e.g. ('C', 'C')
        if pair[0] == pair[1]:
            i1 = np.array(i1)
            j1 = np.array(j1)
            # wrong, offset1 could be non zero
            mask = np.where(
                (array1['indices'][i1] > array2['indices'][j1]), False, True)
            i1 = i1[mask]
            j1 = j1[mask]
            i1 = list(i1)
            j1 = list(j1)
        n = len(i1)
        k1 = [setting[pair]['search']]*n
        p1 = [setting[pair]['polyhedra']]*n
        t1 = [setting[pair]['type']]*n
        i.extend(i1)
        j.extend(j1)
        k.extend(k1)
        p.extend(p1)
        t.extend(t1)
    # offsets
    offsets_i = array1['offsets'][i]
    offsets_j = array2['offsets'][j]
    distance_vector = array1['positions'][i] - array2['positions'][j]
    distances = np.sqrt(np.sum(distance_vector*distance_vector, axis=1))
    #
    i = array1['indices'][i]
    j = array2['indices'][j]
    k = np.array(k)
    p = np.array(p)
    t = np.array(t)
    # Remove all self-interaction.
    if not self_interaction:
        mask = np.where((i == j) &
                        ((array1['offsets'][i] ==
                         array2['offsets'][j]).all(axis=1)),
                        False, True)
        # print(mask)
        i = i[mask]
        j = j[mask]
        k = k[mask]
        p = p[mask]
        t = t[mask]
        distances = distances[mask]
        offsets_i = offsets_i[mask]
        offsets_j = offsets_j[mask]
    #
    retvals = []
    for q in quantities:
        if q == 'i':
            retvals += [i]
        elif q == 'j':
            retvals += [j]
        elif q == 'k':
            retvals += [k]
        elif q == 'p':
            retvals += [p]
        elif q == 't':
            retvals += [t]
        elif q == 'd':
            retvals += [distances]
        elif q == 'S':
            retvals += [offsets_j]
        else:
            raise ValueError('Unsupported quantity specified.')
    if len(retvals) == 1:
        return retvals[0]
    else:
        return tuple(retvals)


def neighbor_kdtree(species0, positions0, cell, pbc,
                    cutoffs):
    """
    wrap to pbc structure

    """
    natom = len(positions0)
    # orignal atoms
    array1 = {
        'positions': positions0,
        'species': species0,
        'indices': np.arange(natom),
        'offsets': np.zeros((natom, 3)),
    }
    # atoms added with boundary
    array2 = RemovePbc(species0, positions0, cell, pbc, cutoffs)
    bonddatas = primitive_neighbor_kdtree(array1, array2, cutoffs)
    # bondlists = neighbor_kdtree(quantities,
    # array1, array2, cutoffs, self_interaction)
    return bonddatas


def primitive_neighbor_kdtree(array1, array2,
                              cutoffs, parallel=1):
    """non pbc
    build bond lists between atoms1 and atoms2.
    In the non-pbc case: atoms1 and atoms2 could be the same.
    In pbc case, atoms2 is atoms1 + boundary atoms
    """
    from scipy.spatial import KDTree
    tstart = time()
    bonddatas = {}
    for pair, cutoff in cutoffs.items():
        bonddatas[pair] = {}
        indices_i = np.where(array1['species'] == pair[0])[0]
        indices_j = np.where(array2['species'] == pair[1])[0]
        if len(indices_i) == 0 or len(indices_j) == 0:
            continue
        indices_min = None
        # qurey the smaller array is faster, flip is needed
        p1 = array1['positions'][indices_i]
        p2 = array2['positions'][indices_j]
        # find max
        tree = KDTree(p2)
        indices_max = tree.query_ball_point(p1, r=cutoff[1], workers=parallel)
        # find min
        if cutoff[0] > 1e-6:
            indices_min = tree.query_ball_point(
                p1, r=cutoff[0], workers=parallel)
        #
        n = len(p1)
        if indices_min is not None:
            for k in range(n):
                indices_max[k] = list(
                    set(indices_max[k]) - set(indices_min[k]))
                m = len(indices_max[k])
                if m == 0:
                    continue
                bonddatas[pair][indices_i[k]] = indices_j[indices_max[k]]
        else:
            for k in range(n):
                m = len(indices_max[k])
                if m == 0:
                    continue
                bonddatas[pair][indices_i[k]] = indices_j[indices_max[k]]
    logger.debug('Build bondlist: {:1.2f}'.format(time() - tstart))
    # print(bonddatas)
    return bonddatas


def cellPlanes(cell, origin=np.array([0, 0, 0])):
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


def pointCellDistance(points, cell, origin=np.zeros(3)):
    """
    get distance between point and the six faces of a cell.
    """
    planes = cellPlanes(cell, origin=origin)
    distances = []
    npoint = len(points)
    distances = np.zeros((3, 2, npoint))
    for i in range(3):
        for j in range(2):
            vec = points - planes[i][j]['points']
            distances[i, j, :] = np.dot(vec, planes[i][j]['normals'])
    return distances


def build_boundary(species, positions,
                   cell, pbc, boundary, include_self=False):
    """
    find region outside cell for periodic caculation, supercell

    todo: boundary > 1
    """
    from functools import reduce
    # tstart = time()
    wraped_positions = wrap_positions(positions, cell, pbc=pbc)
    distances = pointCellDistance(wraped_positions, cell)
    natom = len(positions)
    # find atoms close to cell face with distance of r
    indices = [[[], []], [[], []], [[], []]]
    nb = 0
    for i in range(3):
        if not pbc[i]:
            continue
        ind0 = np.where((distances[i, 0, :] < boundary[i, 1]))[0]
        ind1 = np.where((distances[i, 1, :] < boundary[i, 0]))[0]
        indices[i] = [ind0, ind1]
        nb += len(ind0) + len(ind1)
    # print(indices)
    # init
    # repeat the scaled_positions so that it completely covers the boundary
    positions_b = np.zeros((27*nb, 3))
    offsets_b = np.zeros((27*nb, 3))
    indices_b = np.zeros(27*nb, dtype=int)
    species_b = np.zeros(27*nb, dtype='U20')
    # build face
    offset = [1, -1]
    nb1 = 0
    nb2 = 0
    for c in range(3):
        if not pbc[i]:
            continue
        for i in range(2):
            nb2 = nb1 + len(indices[c][i])
            indices_b[nb1:nb2] = indices[c][i]
            positions_b[nb1:nb2] = positions[indices[c][i]]
            offsets_b[nb1:nb2][:, c] = offset[i]
            nb1 = nb2
    # build edge
    for c in [[0, 1], [0, 2], [1, 2]]:
        if not pbc[c[0]] or not pbc[c[1]]:
            continue
        for i in range(2):
            for j in range(2):
                indices_j = np.intersect1d(indices[c[0]][i], indices[c[1]][j])
                nb2 = nb1 + len(indices_j)
                indices_b[nb1:nb2] = indices_j
                positions_b[nb1:nb2] = positions[indices_j]
                offsets_b[nb1:nb2][:, c[0]] = offset[i]
                offsets_b[nb1:nb2][:, c[1]] = offset[j]
                nb1 = nb2
    # build corner
    if pbc[0] and pbc[1] and pbc[2]:
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    indices3 = reduce(np.intersect1d,
                                      (indices[0][i],
                                       indices[1][j],
                                       indices[2][k]))
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
        positions_b = np.append(positions, positions_b, axis=0)
        indices_b = np.append(np.arange(natom), indices_b)
        species_b = np.append(species, species_b)
        offsets_b = np.append(np.zeros((natom, 3)), offsets_b, axis=0)

    # print('build boundary: {:1.2f}'.format(time() - tstart))
    boundary_data = {
        'positions': positions_b,
        'indices': indices_b,
        'species': species_b,
        'offsets': offsets_b,
    }
    return boundary_data
