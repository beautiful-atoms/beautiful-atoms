import numpy as np
import math
from time import time
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


def map_volumetric_data(volumetric_data, coordinates):
    """Interpolate value at given coordinates for
    the volumetric data

    Note: both volumetric data and coordinate are using scaled positions
    inside the cell
    Args:
        volumetric_data (array): volumetric data
        coordinates (array): coordinates

    Returns:
        array: value at given coordinate
    """
    from scipy import ndimage
    # get scaled coordinates
    index = coordinates*volumetric_data.shape
    # map new value
    data = ndimage.map_coordinates(volumetric_data, index.T, order=1)
    return data

def map_color(data, color1 = [1, 0, 0, 1], color2=[0, 0, 1, 1]):
    """_summary_

    Args:
        data (_type_): _description_
        color1 (list, optional): _description_. Defaults to [1, 0, 0, 1].
        color2 (list, optional): _description_. Defaults to [0, 0, 1, 1].

    Returns:
        _type_: _description_
    """
    # normalize
    data = (data - np.min(data))/(np.max(data) - np.min(data))
    # generate color based on value nvc
    color1 = np.array(color1)
    color2 = np.array(color2)
    dcolor = (color2 - color1)
    color_array = data[:, None]*dcolor
    color = color1 + color_array
    return color

def read_from_others(from_ase=None, from_pymatgen=None,
                     from_pybel=None):
    if from_ase is not None:
        return read_from_ase(from_ase)
    elif from_pymatgen is not None:
        return read_from_pymatgen(from_pymatgen)
    elif from_pybel is not None:
        return read_from_pybel(from_pybel)
    else:
        return None


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
    arrays = atoms.arrays
    if 'species' not in arrays:
        species = np.array(atoms.get_chemical_symbols(), dtype='U20')
    else:
        species = arrays['species']
    info = atoms.info
    # read frames
    if nframe > 1:
        positions = np.zeros((nframe, natom, 3))
        for i in range(0, nframe):
            positions[i, :, :] = frames[i].positions
    else:
        positions = arrays['positions']
    # attributes
    attributes = {}
    for key, array in arrays.items():
        if key in ["positions", "species"]:
            continue
        attributes[key] = array=array
    return species, positions, attributes, atoms.cell, \
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


def read_from_pybel(mol):
    """
    conda install openbabel -c conda-forge
    Converts a pybel molecule to batoms.
    Molecules have the following attributes: atoms, charge, data, dim, energy, exactmass, formula, molwt, spin, sssr, title and unitcell (if crystal data).

    Atoms have the following attributes: atomicmass, atomicnum, coords, exactmass, formalcharge, heavyvalence, heterovalence, hyb, idx, implicitvalence, isotope, partialcharge, spin, type, valence, vector. The .coords attribute provides a tuple (x, y, z) of the atom’s coordinates. The remaining attributes are as for the Get methods of OBAtom.
    """
    from openbabel import pybel
    from ase.data import chemical_symbols
    number = []
    positions = []
    charges = []
    elements = []
    for atom in mol.atoms:
        number.append(atom.atomicnum)
        elements.append(chemical_symbols[atom.atomicnum])
        positions.append(atom.coords)
        charges.append(atom.formalcharge)
    cell = mol.unitcell if hasattr(mol, "unitcell") else [0, 0, 0]
    pbc = True if hasattr(mol, "unitcell") else False
    species = elements.copy()
    arrays = {'charges': charges}
    bonds = [{'atoms': [b.GetBeginAtom().GetIndex(),
                        b.GetEndAtom().GetIndex()],
              'order': b.GetBondOrder()}
             for b in pybel.ob.OBMolBondIter(mol.OBMol)]
    return species, positions, arrays, cell, pbc, {'bonds': bonds}


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
    """Set default color, radii for elements,

    Args:
        elements (dict): _description_
        radius_style (str, optional): _description_. Defaults to 'covalent'.
        color_style (str, optional): _description_. Defaults to "JMOL".
        props (dict, optional): _description_. Defaults to {}.

    exmaples:
    elements = {"Fe": 0.8, "Cr": 0.2}
    elements = {"Fe": {"occupancy": 0.8}, "Cr": {"occupancy": 0.2}}

    Returns:
        _type_: _description_
    """
    radius = 0
    species_props = {"elements": {}}
    for ele, eledata in elements.items():
        if isinstance(eledata, (int, float)):
            # only has occupancy data
            eledata = {"occupancy": eledata}
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

def getDistances(points1, points2):
    """Calculate the distance between every point pair in points.
    Parameters
    ----------
    points1: array_like, shape (n1, m)
        m is the dimension of point.
    points2: array_like, shape (n2, m)
        m is the dimension of point.
    Returns
    -------
    indices : array of integers, shape (n1*n2, 2)
        The index of each point in the pairs.
    dist : array of floats, shape (n1*n2)
        The distances of the point pairs.
    """
    # tstart = time()
    npoint1 = len(points1)
    npoint2 = len(points2)
    i, j = np.triu_indices(npoint1, 0, npoint2)
    vectors = points2[j] - points1[i]
    dist = np.linalg.norm(vectors, axis = 1)
    indices = np.column_stack([i, j])
    # print("getDistance: {}".format(time() - tstart))
    return indices, dist

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




def get_equivalent_indices(no, indices):
    """
    """
    from ase.spacegroup import Spacegroup
    sg = Spacegroup(no)
    indices = sg.equivalent_reflections([indices])
    indices = indices.tolist()
    return indices

def get_equivalent_atoms(atoms, tol = 1e-5):
    """Get equivalent atoms using spglib.
    """
    import spglib
    atoms = (atoms.get_cell(), atoms.get_scaled_positions(),
                   atoms.numbers)
    symmetry_data = spglib.get_symmetry_dataset(atoms, symprec=tol)
    return symmetry_data['equivalent_atoms']

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
    logger.debug('calc_euler_angle: {0:10.2f} s'.format(time() - tstart))
    return eulers

def type_blender_to_py(dtype, str = "str"):
    """Change Blender data type to python data type

    Args:
        type (_type_): _description_
    """
    type_dict = {
        "INT":"int",
        "FLOAT":"float",
        "STRING":str,
        "BOOLEAN":"bool",
    }
    return type_dict[dtype]

def type_py_to_blender(dtype):
    """change python data type to Blender data type

    Args:
        dtype (_type_): _description_

    Returns:
        _type_: _description_
    """
    # Update the subdtype check according to
        # https://github.com/numpy/numpy/blob/db481babcfa7ebc70833e77985858e9295a3135b/numpy/core/numerictypes.py#L357
    if np.issubdtype(dtype, np.integer):
        dtype = 'INT'
    elif np.issubdtype(dtype, np.floating):
        dtype = 'FLOAT'
    elif np.issubdtype(dtype, np.str_):
        dtype = 'STRING'
    elif np.issubdtype(dtype, np.bool_):
        dtype = 'BOOLEAN'
    else:
        # raise KeyError('%s is not supported.' % dtype)
        # print('Attribute: {}, {} is not supported.'.format(dtype))
        return False
    return dtype

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

def deprecated(msg):
    """Helper function print deprecating warning.
    """
    logger.warning(msg)
    msg = "="*80 + "\n" + "Warning: " + msg + "\n" + "="*80 + "\n"
    print(msg)
