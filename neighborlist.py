"""
From ASE
https://wiki.fysik.dtu.dk/ase/_modules/ase/neighborlist.html#neighbor_list

Modified to support:

1) species instead of element
2) add a bondlengt range
    e.g cutoff = {('Ti_1', 'O'):[0.5, 3.0] }

"""

import numpy as np

from ase.data import atomic_numbers
from ase.geometry import complete_cell


import time

def primitive_neighbor_list(quantities, pbc, cell, positions, cutoff,
                            species=None, self_interaction=False,
                            use_scaled_positions=False, max_nbins=1e6):
    """Compute a neighbor list for an atomic configuration.

    Atoms outside periodic boundaries are mapped into the box. Atoms
    outside nonperiodic boundaries are included in the neighbor list
    but complexity of neighbor list search for those can become n^2.

    The neighbor list is sorted by first atom index 'i', but not by second
    atom index 'j'.

    Parameters:

    quantities: str
        Quantities to compute by the neighbor list algorithm. Each character
        in this string defines a quantity. They are returned in a tuple of
        the same order. Possible quantities are

            * 'i' : first atom index
            * 'j' : second atom index
            * 'd' : absolute distance
            * 'D' : distance vector
            * 'S' : shift vector (number of cell boundaries crossed by the bond
              between atom i and j). With the shift vector S, the
              distances D between atoms can be computed from:
              D = positions[j]-positions[i]+S.dot(cell)
    pbc: array_like
        3-tuple indicating giving periodic boundaries in the three Cartesian
        directions.
    cell: 3x3 matrix
        Unit cell vectors.
    positions: list of xyz-positions
        Atomic positions.  Anything that can be converted to an ndarray of
        shape (n, 3) will do: [(x1,y1,z1), (x2,y2,z2), ...]. If
        use_scaled_positions is set to true, this must be scaled positions.
    cutoff: float or dict
        Cutoff for neighbor search. It can be:

            * A single float: This is a global cutoff for all elements.
            * A dictionary: This specifies cutoff values for element
              pairs. Specification accepts element numbers of symbols.
              Example: {(1, 6): 1.1, (1, 1): 1.0, ('C', 'C'): 1.85}
            * A list/array with a per atom value: This specifies the radius of
              an atomic sphere for each atoms. If spheres overlap, atoms are
              within each others neighborhood. See :func:`~ase.neighborlist.natural_cutoffs`
              for an example on how to get such a list.
    self_interaction: bool
        Return the atom itself as its own neighbor if set to true.
        Default: False
    use_scaled_positions: bool
        If set to true, positions are expected to be scaled positions.
    max_nbins: int
        Maximum number of bins used in neighbor search. This is used to limit
        the maximum amount of memory required by the neighbor list.

    Returns:

    i, j, ... : array
        Tuple with arrays for each quantity specified above. Indices in `i`
        are returned in ascending order 0..len(a)-1, but the order of (i,j)
        pairs is not guaranteed.

    """

    # Naming conventions: Suffixes indicate the dimension of an array. The
    # following convention is used here:
    #     c: Cartesian index, can have values 0, 1, 2
    #     i: Global atom index, can have values 0..len(a)-1
    #     xyz: Bin index, three values identifying x-, y- and z-component of a
    #          spatial bin that is used to make neighbor search O(n)
    #     b: Linearized version of the 'xyz' bin index
    #     a: Bin-local atom index, i.e. index identifying an atom *within* a
    #        bin
    #     p: Pair index, can have value 0 or 1
    #     n: (Linear) neighbor index

    # Return empty neighbor list if no atoms are passed here
    tstart = time.time()
    if len(positions) == 0:
        empty_types = dict(i=(int, (0, )),
                           j=(int, (0, )),
                           D=(float, (0, 3)),
                           d=(float, (0, )),
                           S=(int, (0, 3)))
        retvals = []
        for i in quantities:
            dtype, shape = empty_types[i]
            retvals += [np.array([], dtype=dtype).reshape(shape)]
        if len(retvals) == 1:
            return retvals[0]
        else:
            return tuple(retvals)

    # Compute reciprocal lattice vectors.
    b1_c, b2_c, b3_c = np.linalg.pinv(cell).T

    # Compute distances of cell faces.
    l1 = np.linalg.norm(b1_c)
    l2 = np.linalg.norm(b2_c)
    l3 = np.linalg.norm(b3_c)
    face_dist_c = np.array([1 / l1 if l1 > 0 else 1,
                            1 / l2 if l2 > 0 else 1,
                            1 / l3 if l3 > 0 else 1])

    if isinstance(cutoff, dict):
        max_cutoff = np.array([x for x in cutoff.values()]).max()
    else:
        if np.isscalar(cutoff):
            max_cutoff = cutoff
        else:
            cutoff = np.asarray(cutoff)
            max_cutoff = 2*np.max(cutoff)

    # We use a minimum bin size of 3 A
    bin_size = max(max_cutoff, 3)
    # Compute number of bins such that a sphere of radius cutoff fits into
    # eight neighboring bins.
    nbins_c = np.maximum((face_dist_c / bin_size).astype(int), [1, 1, 1])
    nbins = np.prod(nbins_c)
    # Make sure we limit the amount of memory used by the explicit bins.
    while nbins > max_nbins:
        nbins_c = np.maximum(nbins_c // 2, [1, 1, 1])
        nbins = np.prod(nbins_c)

    # Compute over how many bins we need to loop in the neighbor list search.
    neigh_search_x, neigh_search_y, neigh_search_z = \
        np.ceil(bin_size * nbins_c / face_dist_c).astype(int)

    # If we only have a single bin and the system is not periodic, then we
    # do not need to search neighboring bins
    neigh_search_x = 0 if nbins_c[0] == 1 and not pbc[0] else neigh_search_x
    neigh_search_y = 0 if nbins_c[1] == 1 and not pbc[1] else neigh_search_y
    neigh_search_z = 0 if nbins_c[2] == 1 and not pbc[2] else neigh_search_z

    # Sort atoms into bins.
    if use_scaled_positions:
        scaled_positions_ic = positions
        positions = np.dot(scaled_positions_ic, cell)
    else:
        scaled_positions_ic = np.linalg.solve(complete_cell(cell).T,
                                              positions.T).T
    bin_index_ic = np.floor(scaled_positions_ic*nbins_c).astype(int)
    cell_shift_ic = np.zeros_like(bin_index_ic)

    for c in range(3):
        if pbc[c]:
            # (Note: np.divmod does not exist in older numpies)
            cell_shift_ic[:, c], bin_index_ic[:, c] = \
                divmod(bin_index_ic[:, c], nbins_c[c])
        else:
            bin_index_ic[:, c] = np.clip(bin_index_ic[:, c], 0, nbins_c[c]-1)

    # Convert Cartesian bin index to unique scalar bin index.
    bin_index_i = (bin_index_ic[:, 0] +
                   nbins_c[0] * (bin_index_ic[:, 1] +
                                 nbins_c[1] * bin_index_ic[:, 2]))

    # atom_i contains atom index in new sort order.
    atom_i = np.argsort(bin_index_i)
    bin_index_i = bin_index_i[atom_i]

    # Find max number of atoms per bin
    max_natoms_per_bin = np.bincount(bin_index_i).max()

    # Sort atoms into bins: atoms_in_bin_ba contains for each bin (identified
    # by its scalar bin index) a list of atoms inside that bin. This list is
    # homogeneous, i.e. has the same size *max_natoms_per_bin* for all bins.
    # The list is padded with -1 values.
    atoms_in_bin_ba = -np.ones([nbins, max_natoms_per_bin], dtype=int)
    for i in range(max_natoms_per_bin):
        # Create a mask array that identifies the first atom of each bin.
        mask = np.append([True], bin_index_i[:-1] != bin_index_i[1:])
        # Assign all first atoms.
        atoms_in_bin_ba[bin_index_i[mask], i] = atom_i[mask]

        # Remove atoms that we just sorted into atoms_in_bin_ba. The next
        # "first" atom will be the second and so on.
        mask = np.logical_not(mask)
        atom_i = atom_i[mask]
        bin_index_i = bin_index_i[mask]

    # Make sure that all atoms have been sorted into bins.
    assert len(atom_i) == 0
    assert len(bin_index_i) == 0

    # Now we construct neighbor pairs by pairing up all atoms within a bin or
    # between bin and neighboring bin. atom_pairs_pn is a helper buffer that
    # contains all potential pairs of atoms between two bins, i.e. it is a list
    # of length max_natoms_per_bin**2.
    atom_pairs_pn = np.indices((max_natoms_per_bin, max_natoms_per_bin),
                               dtype=int)
    atom_pairs_pn = atom_pairs_pn.reshape(2, -1)

    # Initialized empty neighbor list buffers.
    first_at_neightuple_nn = []
    secnd_at_neightuple_nn = []
    cell_shift_vector_x_n = []
    cell_shift_vector_y_n = []
    cell_shift_vector_z_n = []

    # This is the main neighbor list search. We loop over neighboring bins and
    # then construct all possible pairs of atoms between two bins, assuming
    # that each bin contains exactly max_natoms_per_bin atoms. We then throw
    # out pairs involving pad atoms with atom index -1 below.
    binz_xyz, biny_xyz, binx_xyz = np.meshgrid(np.arange(nbins_c[2]),
                                               np.arange(nbins_c[1]),
                                               np.arange(nbins_c[0]),
                                               indexing='ij')
    # The memory layout of binx_xyz, biny_xyz, binz_xyz is such that computing
    # the respective bin index leads to a linearly increasing consecutive list.
    # The following assert statement succeeds:
    #     b_b = (binx_xyz + nbins_c[0] * (biny_xyz + nbins_c[1] *
    #                                     binz_xyz)).ravel()
    #     assert (b_b == np.arange(np.prod(nbins_c))).all()


    # print('neighborlist 0: {0:10.2f} s'.format( time.time() - tstart))
    tstart = time.time()

    # First atoms in pair.
    _first_at_neightuple_n = atoms_in_bin_ba[:, atom_pairs_pn[0]]
    for dz in range(-neigh_search_z, neigh_search_z+1):
        for dy in range(-neigh_search_y, neigh_search_y+1):
            for dx in range(-neigh_search_x, neigh_search_x+1):
                # Bin index of neighboring bin and shift vector.
                shiftx_xyz, neighbinx_xyz = divmod(binx_xyz + dx, nbins_c[0])
                shifty_xyz, neighbiny_xyz = divmod(biny_xyz + dy, nbins_c[1])
                shiftz_xyz, neighbinz_xyz = divmod(binz_xyz + dz, nbins_c[2])
                neighbin_b = (neighbinx_xyz + nbins_c[0] *
                              (neighbiny_xyz + nbins_c[1] * neighbinz_xyz)
                              ).ravel()

                # Second atom in pair.
                _secnd_at_neightuple_n = \
                    atoms_in_bin_ba[neighbin_b][:, atom_pairs_pn[1]]

                # Shift vectors.
                _cell_shift_vector_x_n = \
                    np.resize(shiftx_xyz.reshape(-1, 1),
                              (max_natoms_per_bin**2, shiftx_xyz.size)).T
                _cell_shift_vector_y_n = \
                    np.resize(shifty_xyz.reshape(-1, 1),
                              (max_natoms_per_bin**2, shifty_xyz.size)).T
                _cell_shift_vector_z_n = \
                    np.resize(shiftz_xyz.reshape(-1, 1),
                              (max_natoms_per_bin**2, shiftz_xyz.size)).T

                # We have created too many pairs because we assumed each bin
                # has exactly max_natoms_per_bin atoms. Remove all surperfluous
                # pairs. Those are pairs that involve an atom with index -1.
                mask = np.logical_and(_first_at_neightuple_n != -1,
                                      _secnd_at_neightuple_n != -1)
                if mask.sum() > 0:
                    first_at_neightuple_nn += [_first_at_neightuple_n[mask]]
                    secnd_at_neightuple_nn += [_secnd_at_neightuple_n[mask]]
                    cell_shift_vector_x_n += [_cell_shift_vector_x_n[mask]]
                    cell_shift_vector_y_n += [_cell_shift_vector_y_n[mask]]
                    cell_shift_vector_z_n += [_cell_shift_vector_z_n[mask]]

    # Flatten overall neighbor list.
    first_at_neightuple_n = np.concatenate(first_at_neightuple_nn)
    secnd_at_neightuple_n = np.concatenate(secnd_at_neightuple_nn)
    cell_shift_vector_n = np.transpose([np.concatenate(cell_shift_vector_x_n),
                                        np.concatenate(cell_shift_vector_y_n),
                                        np.concatenate(cell_shift_vector_z_n)])

    # Add global cell shift to shift vectors
    cell_shift_vector_n += cell_shift_ic[first_at_neightuple_n] - \
        cell_shift_ic[secnd_at_neightuple_n]

    # Remove all self-pairs that do not cross the cell boundary.
    if not self_interaction:
        m = np.logical_not(np.logical_and(
            first_at_neightuple_n == secnd_at_neightuple_n,
            (cell_shift_vector_n == 0).all(axis=1)))
        first_at_neightuple_n = first_at_neightuple_n[m]
        secnd_at_neightuple_n = secnd_at_neightuple_n[m]
        cell_shift_vector_n = cell_shift_vector_n[m]

    # For nonperiodic directions, remove any bonds that cross the domain
    # boundary.
    for c in range(3):
        if not pbc[c]:
            m = cell_shift_vector_n[:, c] == 0
            first_at_neightuple_n = first_at_neightuple_n[m]
            secnd_at_neightuple_n = secnd_at_neightuple_n[m]
            cell_shift_vector_n = cell_shift_vector_n[m]

    # Sort neighbor list.
    i = np.argsort(first_at_neightuple_n)
    first_at_neightuple_n = first_at_neightuple_n[i]
    secnd_at_neightuple_n = secnd_at_neightuple_n[i]
    cell_shift_vector_n = cell_shift_vector_n[i]

    # Compute distance vectors.
    distance_vector_nc = positions[secnd_at_neightuple_n] - \
        positions[first_at_neightuple_n] + \
        cell_shift_vector_n.dot(cell)
    abs_distance_vector_n = \
        np.sqrt(np.sum(distance_vector_nc*distance_vector_nc, axis=1))

    # We have still created too many pairs. Only keep those with distance
    # smaller than max_cutoff.
    mask = abs_distance_vector_n < max_cutoff
    first_at_neightuple_n = first_at_neightuple_n[mask]
    secnd_at_neightuple_n = secnd_at_neightuple_n[mask]
    cell_shift_vector_n = cell_shift_vector_n[mask]
    distance_vector_nc = distance_vector_nc[mask]
    abs_distance_vector_n = abs_distance_vector_n[mask]
    # print('neighborlist 1: {0:10.2f} s'.format( time.time() - tstart))

    tstart = time.time()
    if isinstance(cutoff, dict) and species is not None:
        # If cutoff is a dictionary, then the cutoff radii are specified per
        # element pair. We now have a list up to maximum cutoff.
        species = np.array(species)
        per_pair_cutoff_n = np.zeros_like(abs_distance_vector_n)
        per_pair_cutoff_n_min = np.zeros_like(abs_distance_vector_n)
        for (sp1, sp2), c in cutoff.items():
            mask = np.logical_and(
                    species[first_at_neightuple_n] == sp1,
                    species[secnd_at_neightuple_n] == sp2)
            per_pair_cutoff_n_min[mask] = c[0]
            per_pair_cutoff_n[mask] = c[1]
        mask1 = abs_distance_vector_n < per_pair_cutoff_n
        mask2 = abs_distance_vector_n > per_pair_cutoff_n_min
        mask = mask1 & mask2
        first_at_neightuple_n = first_at_neightuple_n[mask]
        secnd_at_neightuple_n = secnd_at_neightuple_n[mask]
        cell_shift_vector_n = cell_shift_vector_n[mask]
        distance_vector_nc = distance_vector_nc[mask]
        abs_distance_vector_n = abs_distance_vector_n[mask]
    elif not np.isscalar(cutoff):
        # If cutoff is neither a dictionary nor a scalar, then we assume it is
        # a list or numpy array that contains atomic radii. Atoms are neighbors
        # if their radii overlap.
        mask = abs_distance_vector_n < \
            cutoff[first_at_neightuple_n] + cutoff[secnd_at_neightuple_n]
        first_at_neightuple_n = first_at_neightuple_n[mask]
        secnd_at_neightuple_n = secnd_at_neightuple_n[mask]
        cell_shift_vector_n = cell_shift_vector_n[mask]
        distance_vector_nc = distance_vector_nc[mask]
        abs_distance_vector_n = abs_distance_vector_n[mask]
    # print('neighborlist 2: {0:10.2f} s'.format( time.time() - tstart))
    # Assemble return tuple.
    retvals = []
    for q in quantities:
        if q == 'i':
            retvals += [first_at_neightuple_n]
        elif q == 'j':
            retvals += [secnd_at_neightuple_n]
        elif q == 'D':
            retvals += [distance_vector_nc]
        elif q == 'd':
            retvals += [abs_distance_vector_n]
        elif q == 'S':
            retvals += [cell_shift_vector_n]
        else:
            raise ValueError('Unsupported quantity specified.')
    if len(retvals) == 1:
        return retvals[0]
    else:
        return tuple(retvals)



def neighbor_list(quantities, a, cutoff, self_interaction=False,
                  max_nbins=1e6):
    """Compute a neighbor list for an atomic configuration.

    Atoms outside periodic boundaries are mapped into the box. Atoms
    outside nonperiodic boundaries are included in the neighbor list
    but complexity of neighbor list search for those can become n^2.

    The neighbor list is sorted by first atom index 'i', but not by second
    atom index 'j'.

    Parameters:

    quantities: str
        Quantities to compute by the neighbor list algorithm. Each character
        in this string defines a quantity. They are returned in a tuple of
        the same order. Possible quantities are:

           * 'i' : first atom index
           * 'j' : second atom index
           * 'd' : absolute distance
           * 'D' : distance vector
           * 'S' : shift vector (number of cell boundaries crossed by the bond
             between atom i and j). With the shift vector S, the
             distances D between atoms can be computed from:
             D = a.positions[j]-a.positions[i]+S.dot(a.cell)
    a: :class:`ase.Atoms`
        Atomic configuration.
    cutoff: float or dict
        Cutoff for neighbor search. It can be:

            * A single float: This is a global cutoff for all elements.
            * A dictionary: This specifies cutoff values for element
              pairs. Specification accepts element numbers of symbols.
              Example: {(1, 6): 1.1, (1, 1): 1.0, ('C', 'C'): 1.85}
            * A list/array with a per atom value: This specifies the radius of
              an atomic sphere for each atoms. If spheres overlap, atoms are
              within each others neighborhood. See :func:`~ase.neighborlist.natural_cutoffs`
              for an example on how to get such a list.

    self_interaction: bool
        Return the atom itself as its own neighbor if set to true.
        Default: False
    max_nbins: int
        Maximum number of bins used in neighbor search. This is used to limit
        the maximum amount of memory required by the neighbor list.

    Returns:

    i, j, ...: array
        Tuple with arrays for each quantity specified above. Indices in `i`
        are returned in ascending order 0..len(a), but the order of (i,j)
        pairs is not guaranteed.

    Examples:

    Examples assume Atoms object *a* and numpy imported as *np*.

    1. Coordination counting::

        i = neighbor_list('i', a, 1.85)
        coord = np.bincount(i)

    2. Coordination counting with different cutoffs for each pair of species::

        i = neighbor_list('i', a,
                          {('H', 'H'): 1.1, ('C', 'H'): 1.3, ('C', 'C'): 1.85})
        coord = np.bincount(i)

    3. Pair distribution function::

        d = neighbor_list('d', a, 10.00)
        h, bin_edges = np.histogram(d, bins=100)
        pdf = h/(4*np.pi/3*(bin_edges[1:]**3 - bin_edges[:-1]**3)) * a.get_volume()/len(a)

    4. Pair potential::

        i, j, d, D = neighbor_list('ijdD', a, 5.0)
        energy = (-C/d**6).sum()
        pair_forces = (6*C/d**5  * (D/d).T).T
        forces_x = np.bincount(j, weights=pair_forces[:, 0], minlength=len(a)) - \
                   np.bincount(i, weights=pair_forces[:, 0], minlength=len(a))
        forces_y = np.bincount(j, weights=pair_forces[:, 1], minlength=len(a)) - \
                   np.bincount(i, weights=pair_forces[:, 1], minlength=len(a))
        forces_z = np.bincount(j, weights=pair_forces[:, 2], minlength=len(a)) - \
                   np.bincount(i, weights=pair_forces[:, 2], minlength=len(a))

    5. Dynamical matrix for a pair potential stored in a block sparse format::

        from scipy.sparse import bsr_matrix
        i, j, dr, abs_dr = neighbor_list('ijDd', atoms)
        energy = (dr.T / abs_dr).T
        dynmat = -(dde * (energy.reshape(-1, 3, 1) * energy.reshape(-1, 1, 3)).T).T \
                 -(de / abs_dr * (np.eye(3, dtype=energy.dtype) - \
                   (energy.reshape(-1, 3, 1) * energy.reshape(-1, 1, 3))).T).T
        dynmat_bsr = bsr_matrix((dynmat, j, first_i), shape=(3*len(a), 3*len(a)))

        dynmat_diag = np.empty((len(a), 3, 3))
        for x in range(3):
            for y in range(3):
                dynmat_diag[:, x, y] = -np.bincount(i, weights=dynmat[:, x, y])

        dynmat_bsr += bsr_matrix((dynmat_diag, np.arange(len(a)),
                                  np.arange(len(a) + 1)),
                                 shape=(3 * len(a), 3 * len(a)))

    """
    return primitive_neighbor_list(quantities, a.pbc,
                                   a.get_cell(complete=True),
                                   a.positions, cutoff, species=a.arrays['species'],
                                   self_interaction=self_interaction,
                                   max_nbins=max_nbins)

def neighbor_kdtree(quantities, species, positions, cell, pbc,
                    cutoffs, self_interaction=False,
                    use_scaled_positions=False,
                    max_nbins=1e6):
    """
    
    """
    tstart1 = time.time()
    if use_scaled_positions:
        scaled_positions = positions
        positions = np.dot(scaled_positions, cell)
    else:
        scaled_positions = np.linalg.solve(complete_cell(cell).T,
                                              positions.T).T
    # find region outside cell for periodic caculation, supercell 
    max_cutoff = np.array([x for x in cutoffs.values()]).max()
    # Compute distances of cell faces.
    face_dist_c = 1/(np.linalg.norm(np.linalg.pinv(cell).T, axis = 1) + 1e-9)
    boundary = max_cutoff/face_dist_c
    ib = np.array([np.floor(-boundary), np.ceil(1+boundary)], dtype=int)
    for c in range(3):
        if not pbc[c]:
            ib[:, c] = [0, 1]
    M = np.product(ib[1] - ib[0])
    natom = len(scaled_positions)
    # repeat the scaled_positions so that it completely covers the boundary
    scaled_positions_b = np.tile(scaled_positions, (M, ) + (1, ) * (len(scaled_positions.shape) - 1))
    offsets_b = np.zeros((M*natom, 3), dtype=int)
    tstart1 = time.time()
    indices0 = np.arange(natom)
    indices_b = np.tile(indices0, (M, ) + (1, ) * (len(indices0.shape) - 1))
    species_b = np.tile(species, (M, ) + (1, ) * (len(species.shape) - 1))
    i0 = 0
    for m0 in range(ib[0, 0], ib[1, 0]):
        for m1 in range(ib[0, 1], ib[1, 1]):
            for m2 in range(ib[0, 2], ib[1, 2]):
                i1 = i0 + natom
                scaled_positions_b[i0:i1] += np.array([m0, m1, m2])
                offsets_b[i0:i1] = np.array([m0, m1, m2])
                i0 = i1
    # find atoms inside boundary
    ind =  np.where((scaled_positions_b[:, 0] > -boundary[0]) & (scaled_positions_b[:, 0] < 1 + boundary[0]) \
                & (scaled_positions_b[:, 1] > -boundary[1]) & (scaled_positions_b[:, 1] < 1 + boundary[1]) \
                & (scaled_positions_b[:, 2] > -boundary[2]) & (scaled_positions_b[:, 2] < 1 + boundary[2]))[0]
    scaled_positions_b = scaled_positions_b[ind]
    indices_b = indices_b[ind]
    species_b = species_b[ind]
    offsets_b = offsets_b[ind]
    print('build boundary: {:1.2f}'.format(time.time() - tstart1))
    #
    i = []
    j = []
    j_b = []
    positions_b = np.dot(scaled_positions_b, cell)
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
    tstart = time.time()
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
    print('Build distances: {:1.2f}'.format(time.time() - tstart))
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
    if len(retvals) == 1:
        return retvals[0]
    else:
        return tuple(retvals)

def primitive_neighbor_kdtree(positions1, positions2, 
                cutoff = 1.0, parallel = 1):
    """
    """
    from scipy.spatial import KDTree
    #
    tstart = time.time()
    tree = KDTree(positions2)
    indices = tree.query_ball_point(positions1, r = cutoff, workers=parallel)
    print('KDTree: {:1.2f}'.format(time.time() - tstart))
    return indices

if __name__ == "__main__":
    from ase.io import read, write
    from ase.build import molecule, bulk
    from ase.visualize import view
    from ase import Atom, Atoms
    mof = read('test/datas/mof-5.cif')
    pro = read('test/datas/1tim.pdb')
    # kras = kras*[2, 2, 2]
    # kras.center(2)
    h2o = molecule('CH3CH2OH')
    h2o.pbc = True
    h2o.center(3)
    # h2o.cell = [4, 7, 5]
    atoms = h2o
    atoms = atoms*[1, 1, 1]
    atoms.write('h2o-20-20-20.xyz')
    print(len(atoms))
    # view(atoms)
    atoms.arrays['species'] = atoms.get_chemical_symbols()
    print(atoms.arrays['species'])
    # atoms = read('docs/source/_static/datas/mof-5.cif')
    cutoff = {
        # ('O', 'H'): [0.0, 1.261], 
            # ('O', 'O'): [0.0, 1.39],
            # ('H', 'H'): [0.0, 1.39],
            # ('C', 'H'): [0.0, 1.39],
            # ('C', 'O'): [0.0, 1.8],
            ('C', 'C'): [0.0, 1.99],
                }
    tstart = time.time()
    i, j, d, S = neighbor_list('ijdS', atoms, cutoff)
    print('time %s'%(time.time() - tstart))
    print(i, j, S)
    cutoffs = {
        # (30, 8) : [0.000, 2.444],
        # (30, 1) : [0.000, 1.989],
        # (8, 8)  : [0.000, 1.716],
        (8, 1)  : [0.000, 1.261],
        # (6, 8)  : [0.000, 1.846],
        # (6, 6)  : [0.000, 1.976],
        # (1, 1)  : [0.000, 1.076],
        # (6, 1)  : [0.000, 1.391],
        (1, 8)  : [1.200, 2.100],
        }
    tstart = time.time()
    scaled_positions = atoms.get_scaled_positions()
    species = atoms.numbers
    i, j, d, offsets = neighbor_kdtree('ijdS', species, 
                atoms.positions, atoms.get_cell(complete=True), atoms.pbc,
            cutoffs, use_scaled_positions=False)
    print('time %s'%(time.time() - tstart))
    print(i, j, d, offsets)
