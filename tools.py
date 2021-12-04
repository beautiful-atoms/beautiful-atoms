import numpy as np
import time

from numpy.core.numeric import indices

    

def default_element_prop(element, radius_style = 'covalent', color_style = "JMOL"):
    """
    Get color, radii for element.
    """
    from batoms.data import covalent_radii, vdw_radii
    from ase.data import chemical_symbols
    from batoms.data import jmol_colors, cpk_colors, vesta_color
    element_prop = {}
    number = chemical_symbols.index(element)
    color_style = color_style.upper()
    radius_style = radius_style.upper()
    if color_style == 'JMOL':
        color = jmol_colors[number]
    elif color_style == 'CPK':
        color = cpk_colors[number]
    elif color_style == 'VESTA':
        color = vesta_color[element]
    if radius_style == 'COVALENT' or radius_style == '0':
        radius = covalent_radii[number]
    elif radius_style == 'VDW' or radius_style == '1':
        radius = vdw_radii[number]
    elif radius_style == 'IONIC' or radius_style == '2':
        raise Exception('Ionic radii is not supported yet!')
    else:
        raise Exception('radius_style: %s is not supported yet!'%radius_style)
    element_prop['radius'] = radius
    element_prop['color'] = [color[0], color[1], color[2], 1.0]
    return element_prop

def get_default_species_data(elements, radius_style = 'covalent', 
                color_style = "JMOL", props = {}):
    """
    Set default color, radii for elements,
    Todo fraction occupancy
    """
    species_data = {'color':{}}
    radius = 0
    for ele, fraction in elements.items():
        data = default_element_prop(ele, radius_style = radius_style, 
                color_style = color_style)
        radius += data['radius']*fraction
        species_data['color'][ele] = data['color']
    species_data['radius'] = radius
    species_data.update(props)
    return species_data


def get_polyhedra_kind(color, width = 0.01, show_edge = True, props = {}):
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


def get_canvas(vertices, direction = [0, 0 ,1], padding = 1):
    """
    Calculate canvas for camera. 
    Project vertices and cell on the plane of view.
    Find the boudanry of the object on the plane.
    """
    
    # plane
    canvas = np.zeros([2, 3])
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


def find_cage_sphere(cell, positions, radius, step = 1.0, 
                    boundary = [[0, 1], [0, 1], [0, 1]]):
    from ase.cell import Cell
    from scipy.spatial import distance

    cell = Cell(cell)
    a, b, c, alpha, beta, gamma = cell.cellpar()
    na, nb, nc = int(a/step),int(b/step), int(c/step)
    x = np.linspace(boundary[0][0], boundary[0][1], na)
    y = np.linspace(boundary[1][0], boundary[1][1], nb)
    z = np.linspace(boundary[2][0], boundary[2][1], nc)
    positions_v = np.vstack(np.meshgrid(x, y, z)).reshape(3,-1).T
    positions_v = np.dot(positions_v, cell)
    dists = distance.cdist(positions_v, positions)
    # dists<3.0
    flag = np.min(dists, axis = 1) >radius
    return positions_v[flag]


def get_equivalent_indices(no, indices):
    """
    """
    from ase.spacegroup import Spacegroup
    sg = Spacegroup(no)
    indices = sg.equivalent_reflections([indices])
    indices = indices.tolist()
    return indices


def local2global(positions, matrix, reversed = False):
    if reversed:
        matrix = np.linalg.inv(matrix)
    n = len(positions)
    # positions (natom, 3) to (natom, 4)
    positions = np.append(positions, np.ones((n, 1)), axis = 1)
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

if __name__ == "__main__":
    indices = (1, 1, 1)
    indices = get_equivalent_indices(225, indices)
    print(indices)