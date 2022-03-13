import numpy as np
import math
from time import time


def read_from_ase(atoms):
    """
    Import structure from ASE atoms.
    """
    if isinstance(atoms, list):
        frames = atoms
        atoms = frames[0]
    else:
        frames = [atoms]
    nframe = len(frames)
    natom = len(atoms)
    if 'species' not in atoms.arrays:
        atoms.new_array('species', np.array(
            atoms.get_chemical_symbols(), dtype='U20'))
    if 'scale' not in atoms.arrays:
        atoms.new_array('scale', np.ones(len(atoms)))
    arrays = atoms.arrays
    species = arrays.pop('species')
    info = atoms.info
    # read frames
    if nframe > 1:
        positions = np.zeros((nframe, natom, 3))
        for i in range(0, nframe):
            positions[i, :, :] = frames[i].positions
        arrays.pop('positions')
    else:
        positions = arrays.pop('positions')
    return species, positions, arrays, atoms.cell, \
        npbool2bool(atoms.pbc), info


def read_from_pymatgen(structure):
    """
    Import structure from Pymatgen structure.
    """
    attributes = {}
    symbols = np.array([str(site.specie.symbol) for site in structure])
    natom = len(symbols)
    attributes['species'] = np.array(symbols, dtype='U20')
    if hasattr(structure, "lattice"):
        cell = structure.lattice.matrix
        pbc = True
    else:
        cell = np.zeros(3)
        pbc = False
    positions = [structure[i].coords for i in range(natom)]
    info = {}
    return symbols, positions, attributes, cell, pbc, info


def string2Number(s):
    s = str(s)
    return int.from_bytes(s.encode(), 'little')


def number2String(n):
    n = int(n)
    return n.to_bytes(math.ceil(n.bit_length() / 8), 'little').decode()


def default_element_prop(element, radius_style='covalent',
                         color_style="JMOL"):
    """
    Get color, radii for element.
    """
    from batoms.data import covalent_radii, vdw_radii
    from ase.data import chemical_symbols
    from batoms.data import jmol_colors, cpk_colors, vesta_color
    element_prop = {}
    number = chemical_symbols.index(element)
    color_style = str(color_style)
    radius_style = str(radius_style)
    color_style = color_style.upper()
    radius_style = radius_style.upper()
    if color_style == 'JMOL' or color_style == '0':
        color = jmol_colors[number]
    elif color_style == 'CPK' or color_style == '1':
        color = cpk_colors[number]
    elif color_style == 'VESTA' or color_style == '2':
        color = vesta_color[element]
    if radius_style == 'COVALENT' or radius_style == '0':
        radius = covalent_radii[number]
    elif radius_style == 'VDW' or radius_style == '1':
        radius = vdw_radii[number]
    elif radius_style == 'IONIC' or radius_style == '2':
        raise Exception('Ionic radii is not supported yet!')
    else:
        raise Exception('radius_style: %s is not supported yet!' %
                        radius_style)
    element_prop['radius'] = radius
    element_prop['color'] = [color[0], color[1], color[2], 1.0]
    return element_prop


def get_default_species_data(elements, radius_style='covalent',
                             color_style="JMOL", props={}):
    """
    Set default color, radii for elements,
    Todo fraction occupancy
    """
    radius = 0
    species_props = {"elements": {}}
    for ele, eledata in elements.items():
        species_props["elements"][ele] = eledata
        data = default_element_prop(ele, radius_style=radius_style,
                                    color_style=color_style)
        radius += data['radius']*eledata["occupancy"]
        species_props["elements"][ele]["color"] = data["color"]
    species_props['radius'] = radius
    species_props.update(props)
    return species_props


def get_polyhedra_kind(color, width=0.01, show_edge=True,
                       props={}):
    """
    Set initial data for polyhedra.

    """
    polyhedra_kind = {'color': color,
                      'vertices': [],
                      'edges': [],
                      'faces': [],
                      'show_edge': show_edge,
                      }
    polyhedra_kind['edge'] = {'lengths': [],
                              'centers': [],
                              'normals': [],
                              'color': (1.0, 1.0, 1.0, 1.0),
                              'width': width,
                              }
    polyhedra_kind.update(props)
    return polyhedra_kind


def getEquidistantPoints(p1, p2, n):
    points = zip(np.linspace(p1[0], p2[0], n+1),
                 np.linspace(p1[1], p2[1], n+1),
                 np.linspace(p1[2], p2[2], n+1))
    return points


def get_cell_vertices(cell):
    """
    """
    cell = np.array(cell)
    cell = cell.reshape(3, 3)
    cell_vertices = np.empty((2, 2, 2, 3))
    for c1 in range(2):
        for c2 in range(2):
            for c3 in range(2):
                cell_vertices[c1, c2, c3] = np.dot([c1, c2, c3],
                                                   cell)
    cell_vertices.shape = (8, 3)
    return cell_vertices


def get_box(vertices, padding=4):
    """
    Find a minmum box contains all atoms
    """
    # find box
    box = np.zeros((3, 2))
    for i in range(3):
        box[i, 0] = vertices[:, i].min() - padding
        box[i, 1] = vertices[:, i].max() + padding
    return box


def build_grid(box, resolution):
    """
    generate gridpoints
    """
    grids = []
    for i in range(3):
        grids.append(np.arange(box[i, 0], box[i, 1], resolution))
    x, y, z = np.meshgrid(grids[0], grids[1], grids[2],
                          indexing='ij')  # , sparse=True)
    shape = x.shape
    meshgrids = np.c_[x.ravel(), y.ravel(), z.ravel()]
    return meshgrids, shape


def get_canvas(vertices, direction=[0, 0, 1], padding=1):
    """
    Calculate canvas for camera.
    Project vertices and cell on the plane of view.
    Find the boudanry of the object on the plane.
    """

    # plane
    canvas = np.zeros((2, 3))
    frame = rotate_frame(direction)
    #
    projz = np.dot(vertices, frame[2])
    projxy = vertices - projz[:, None]*frame[2]
    projx = np.dot(projxy, frame[0])
    projy = np.dot(projxy, frame[1])
    #
    # find canvas box
    canvas = np.zeros([2, 3])
    canvas[0, 0] = projx.min()
    canvas[1, 0] = projx.max()
    canvas[0, 1] = projy.min()
    canvas[1, 1] = projy.max()
    canvas[0, 2] = projz.min()
    canvas[1, 2] = projz.max()
    canvas[0] -= padding
    canvas[1] += padding
    return canvas


def rotate_frame(direction):
    """
    rotate frame by algin z to direction
    """
    direction = np.array(direction)
    nz = direction/np.linalg.norm(direction)
    nx = np.cross([0, 0, 1], nz) + np.array([1e-6, 0, 0])
    nx = nx/np.linalg.norm(nx)
    ny = np.cross(nz, nx) + np.array([0, 1e-6, 0])
    ny = ny/np.linalg.norm(ny)
    return np.array([nx, ny, nz])


def find_cage_sphere(cell, positions, radius, step=1.0,
                     boundary=[[0, 1], [0, 1], [0, 1]]):
    from ase.cell import Cell
    from scipy.spatial import distance

    cell = Cell(cell)
    a, b, c, alpha, beta, gamma = cell.cellpar()
    na, nb, nc = int(a/step), int(b/step), int(c/step)
    x = np.linspace(boundary[0][0], boundary[0][1], na)
    y = np.linspace(boundary[1][0], boundary[1][1], nb)
    z = np.linspace(boundary[2][0], boundary[2][1], nc)
    positions_v = np.vstack(np.meshgrid(x, y, z)).reshape(3, -1).T
    positions_v = np.dot(positions_v, cell)
    dists = distance.cdist(positions_v, positions)
    # dists<3.0
    flag = np.min(dists, axis=1) > radius
    return positions_v[flag]


def get_equivalent_indices(no, indices):
    """
    """
    from ase.spacegroup import Spacegroup
    sg = Spacegroup(no)
    indices = sg.equivalent_reflections([indices])
    indices = indices.tolist()
    return indices


def local2global(positions, matrix, reversed=False):
    if reversed:
        matrix = np.linalg.inv(matrix)
    n = len(positions)
    # positions (natom, 3) to (natom, 4)
    positions = np.append(positions, np.ones((n, 1)), axis=1)
    # inverse of transformation matrix
    positions = matrix.dot(positions.T).T
    # (natom, 4) back to (natom, 3)
    positions = positions[:, :3]
    return positions


def npbool2bool(pbc):
    """
    """
    newpbc = []
    for i in range(3):
        if pbc[i]:
            newpbc.append(True)
        else:
            newpbc.append(False)
    return newpbc


def heron3(a, b, c):
    s = (a + b + c)/2
    A = np.sqrt(s*(s - a)*(s - b)*(s - c))
    return A


def heron4(u, v, w, U, V, W):
    X = (w - U + v)*(U + v + w)
    x = (U - v + w)*(v - w + U)
    Y = (u - V + w)*(V + w + u)
    y = (V - w + u)*(w - u + V)
    Z = (v - W + u)*(W + u + v)
    z = (W - u + v)*(u - v + W)
    a = np.sqrt(x*Y*Z)
    b = np.sqrt(y*Z*X)
    c = np.sqrt(z*X*Y)
    d = np.sqrt(x*y*z)
    V = np.sqrt((-a + b + c + d)*(a - b + c + d)
                * (a + b - c + d)*(a + b + c - d))
    V = V/192/u/v/w
    return V


def calc_V2(u, v, w, U, V, W):
    a1 = np.square(u)
    a2 = np.square(v)
    a3 = np.square(w)
    a4 = np.square(W)
    a5 = np.square(U)
    a6 = np.square(V)
    V2 = (a1*a5*(a2 + a3 + a4 + a6 - a1 - a5) +
          a2*a6*(a1 + a3 + a4 + a5 - a2 - a6) +
          a3*a4*(a1 + a2 + a5 + a6 - a3 - a4) -
          a1*a2*a4 - a2*a3*a5 - a1*a3*a6 - a4*a5*a6)/144
    return V2


def heron42(u, v, w, U, V, W):
    V2 = calc_V2(u, v, w, U, V, W)
    V2 = np.where(V2 > 0, V2, 1e-6)
    Volume = np.sqrt(V2)
    return Volume


def calc_origin_2(p, p0, p1, r0, r1, r):
    xaxis = p1 - p0
    # xaxis = np.sign(xaxis[:, 0])[:, None]*xaxis
    l0 = r0 + r
    l1 = r1 + r
    l2 = np.linalg.norm(xaxis, axis=1)
    #
    vec = p - p0
    xaxis = xaxis/l2[:, None]
    zaxis = np.cross(xaxis, vec)
    zaxis = zaxis/np.linalg.norm(zaxis, axis=1)[:, None]
    yaxis = np.cross(zaxis, xaxis)
    v = (np.square(l0) + np.square(l2) - np.square(l1))/2/l0/l2
    v[v > 1] = 1
    # print(min(v), max(v))
    angles = np.arccos(v)
    origin = p0 + l0[:, None]*(np.cos(angles)[:, None]
                               * xaxis + np.sin(angles)[:, None]*yaxis)
    return origin


def check_origin_2(p, p0, p1, r0, r1, r, eps):
    origins = calc_origin_2(p, p0, p1, r0, r1, r)
    vec = p - origins
    norms = np.linalg.norm(vec, axis=1)
    indices = np.where(norms < eps)[0]
    return indices, origins[indices]


def check_origin_3(p, p0, p1, p2, r0, r1, r2, r, eps):
    origins = calc_origin_3(p, p0, p1, p2, r0, r1, r2, r)
    vec = p - origins
    norms = np.linalg.norm(vec, axis=1)
    indices = np.where(norms < eps)[0]
    return indices, origins[indices]


def calc_origin_3(p, p0, p1, p2, r0, r1, r2, r):
    l0 = r0 + r
    l1 = r1 + r
    l2 = r2 + r
    l3 = p1 - p0
    l4 = p2 - p0
    l5 = p2 - p1
    l3 = np.linalg.norm(l3, axis=1)
    l4 = np.linalg.norm(l4, axis=1)
    l5 = np.linalg.norm(l5, axis=1)
    xaxis = p1 - p0
    xaxis = xaxis/np.linalg.norm(xaxis, axis=1)[:, None]
    vec1 = p2 - p0
    zaxis = np.cross(xaxis, vec1)
    zaxis = zaxis/np.linalg.norm(zaxis, axis=1)[:, None]
    yaxis = np.cross(zaxis, xaxis)
    # -------------------------
    V = heron42(l0, l1, l2, l5, l4, l3)
    A = heron3(l3, l4, l5)
    h = 3*V/A
    angles1 = np.arccos(
        (np.square(l3) + np.square(l0) - np.square(l1))/2/l3/l0)
    h1 = l0*np.sin(angles1)
    x = l0*np.cos(angles1)
    y = np.sqrt(h1*h1 - h*h)
    # sign of y
    angles2 = np.arccos(
        (np.square(l4) + np.square(l3) - np.square(l5))/2/l4/l3)
    h2 = np.sqrt(x*x + l4*l4 - 2*x*l4*np.cos(angles2))
    h3 = np.sqrt(l2*l2 - h*h)
    projy = np.where(h2 > h3, 1, -1)*y
    # sign of z
    vec = p - p0
    projz = np.sign(np.einsum('ij,ij->i', vec, zaxis))*h
    origin = p0 + x[:, None]*xaxis + \
        projy[:, None]*yaxis + projz[:, None]*zaxis
    return origin


def calc_euler_angle(x, z, seq='xyz'):
    """
    """
    from scipy.spatial.transform import Rotation as R
    tstart = time()
    if len(x.shape) == 1:
        x = x.reshape(-1, 3)
    if len(z.shape) == 1:
        z = z.reshape(-1, 3)
    y = np.cross(z, x, axis=1)
    n = len(z)
    if n == 0:
        eulers = np.zeros((0, 3))
    else:
        x = x.reshape(n, 1, 3)
        y = y.reshape(n, 1, 3)
        z = z.reshape(n, 1, 3)
        mat = np.concatenate((x, y, z), axis=1)
        mat = np.linalg.inv(mat)
        r = R.from_matrix(mat)
        eulers = r.as_euler(seq=seq)  # , degrees = True)
    print('calc_euler_angle: {0:10.2f} s'.format(time() - tstart))
    return eulers


if __name__ == '__main__':
    # V = heron4(6, 7, 8, 9, 10, 11)
    # V3 = heron4(2.59, 2.92, 2.59, 0.96, 1.52, 0.96)
    # print(V, V2, V3)
    p0 = np.array([[-2.5, 0, 0], [0, 0.7632, -0.4770]])
    p1 = np.array([[0, -0.7632, -0.4770], [0, 0, 0.11926]])
    p2 = np.array([[0, 0, 0.11926], [0, -0.7632, -0.4770]])
    r0 = np.array([1.2, 1.2])
    r1 = np.array([1.2, 1.52])
    r2 = np.array([1.52, 1.2])
    r = 1.4
    p = np.array([[0.89120579, 0.98369867, -0.53714401],
                  [-0.89120579, -0.98369867, 0.53714401]])
    p = np.array([[-0.63896561,  0.47235048, -1.45207691],
                  [-0.80025041,  0.32016024, -1.3009696]])
    # origin = calc_origin(p, p0, p1, r0, r1, r)
    origin = calc_origin_3(p, p0, p1, p2, r0, r1, r2, r)
    print(origin)
