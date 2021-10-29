import numpy as np
import time

from numpy.core.numeric import indices

    
def default_element_prop(element, radii_style = 'covalent', color_style = "JMOL"):
    """
    Get color, radii for element.
    """
    from batoms.data import covalent_radii, vdw_radii
    from ase.data import chemical_symbols
    from batoms.data import jmol_colors, cpk_colors, vesta_color
    element_prop = {}
    number = chemical_symbols.index(element)
    if color_style.upper() == 'JMOL':
        color = jmol_colors[number]
    elif color_style.upper() == 'CPK':
        color = cpk_colors[number]
    elif color_style.upper() == 'VESTA':
        color = vesta_color[element]
    if radii_style.upper() == 'COVALENT':
        radius = covalent_radii[number]
    elif radii_style.upper() == 'VDW':
        radius = vdw_radii[number]
    elif radii_style.upper() == 'IONIC':
        raise Exception('Ionic radii is not supported yet!')
    element_prop['radius'] = radius
    element_prop['element'] = element
    element_prop['color'] = [color[0], color[1], color[2], 1.0]
    return element_prop

def get_default_species_data(element, positions = [], radii_style = 'covalent', 
                color_style = "JMOL", scale = [1, 1, 1], props = {}):
    """
    Set default color, radii for element,
    Todo fraction occupancy
    """
    if isinstance(element, str):
        species_data = default_element_prop(element, radii_style = radii_style, 
                    color_style = color_style)
        species_data['scale'] = scale
        species_data['positions'] = positions
        species_data['balltype'] = None
        species_data.update(props)
    elif isinstance(element, dict):
        species_data = {}
        for ele, fraction in element.items():
            species_data[ele] = default_element_prop(ele, radii_style = radii_style, 
                    color_style = color_style)
            species_data[ele]['scale'] = scale
            species_data[ele]['positions'] = positions
            species_data[ele]['balltype'] = None
            species_data[ele].update(props)
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
def get_canvas(positions, cell_vertices, direction = [0, 0 ,1], margin = 1, 
                show_unit_cell = True):
    """
    Calculate canvas for camera. 
    Project positions and cell on the plane of view.
    Find the boudanry of the object on the plane.
    """
    from scipy.spatial.transform import Rotation as R
    canvas = np.zeros([2, 3])
    for i in range(3):
        canvas[0, i] = positions[:, i].min()
        canvas[1, i] = positions[:, i].max()
        if show_unit_cell:
            canvas[0, i] = min(canvas[0, i], cell_vertices[:, i].min())
            canvas[1, i] = max(canvas[1, i], cell_vertices[:, i].max())
    #
    canvas1 = np.zeros([2, 3])
    direction = np.array(direction)
    nz = direction/np.linalg.norm(direction)
    nx = np.cross([0, 0, 1], nz) + np.array([0.000001, 0, 0])
    nx = nx/np.linalg.norm(nx)
    ny = np.cross(nz, nx) + np.array([0, 0.000001, 0])
    ny = ny/np.linalg.norm(ny)
    #
    projxy = positions.copy()
    projx = np.zeros((len(positions), 1))
    projy = np.zeros((len(positions), 1))
    projz = np.dot(positions, nz)
    for i in range(len(positions)):
        projxy[i] = positions[i] - projz[i]*nz
        projx[i] = np.dot(projxy[i], nx)
        projy[i] = np.dot(projxy[i], ny)
    #
    projcellxy = cell_vertices.copy()
    projcellx = np.zeros((len(cell_vertices), 1))
    projcelly = np.zeros((len(cell_vertices), 1))
    for i in range(len(cell_vertices)):
        projcellxy[i] = cell_vertices[i] - np.dot(cell_vertices[i], nz)*nz
        projcellx[i] = np.dot(projcellxy[i], nx)
        projcelly[i] = np.dot(projcellxy[i], ny)
    canvas1[0, 0] = projx.min()
    canvas1[1, 0] = projx.max()
    canvas1[0, 1] = projy.min()
    canvas1[1, 1] = projy.max()
    if show_unit_cell:
        canvas1[0, 0] = min(canvas1[0, 0], projcellx.min())
        canvas1[1, 0] = max(canvas1[1, 0], projcellx.max())
        canvas1[0, 1] = min(canvas1[0, 1], projcelly.min())
        canvas1[1, 1] = max(canvas1[1, 1], projcelly.max())
    canvas1[0, :] -= margin
    canvas1[1, :] += margin
    canvas1[0, 2] = 0
    canvas1[1, 2] = projz.max()
    return canvas, canvas1
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
    
if __name__ == "__main__":
    indices = (1, 1, 1)
    indices = get_equivalent_indices(225, indices)
    print(indices)