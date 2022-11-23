"""Definition of the CrystalShape object in the Batoms package.

"""

import bpy
from time import time
import numpy as np
from batoms.plugins.base import PluginObject
from .setting import CrystalShapeSettings
from batoms.draw import draw_cylinder, draw_surface_from_vertices
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)

class CrystalShape(PluginObject):
    def __init__(self,
                 label=None,
                 location=np.array([0, 0, 0]),
                 batoms=None,
                 ):
        """CrystalShape Class

        Args:
            label (_type_, optional): _description_. Defaults to None.
            location (_type_, optional): _description_. Defaults to np.array([0, 0, 0]).
            batoms (_type_, optional): _description_. Defaults to None.
        """
        #
        self.batoms = batoms
        self.label = label
        PluginObject.__init__(self, label, 'crystalshape', batoms)
        self.settings = CrystalShapeSettings(
            self.label,  parent=self)
        self.settings.bpy_data.active = True

    def build_materials(self, name, color, node_inputs=None,
                        material_style='default'):
        """
        """
        from batoms.material import create_material
        if name in bpy.data.materials:
            mat = bpy.data.materials.get(name)
            bpy.data.materials.remove(mat, do_unlink=True)
        mat = create_material(name,
                              color=color,
                              node_inputs=node_inputs,
                              material_style=material_style,
                              backface_culling=False)
        return mat

    def build_crystal(self, bcell, origin=[0, 0, 0]):
        """
        Build crystal

        no: Int
            spacegroup number
        """
        self.settings.get_symmetry_indices()
        planes = {}
        for p in self.settings:
            normal = np.dot(list(p.indices), bcell.reciprocal)
            normal = normal/np.linalg.norm(normal)
            point = p.distance*normal
            planes[p.name] = {
                'indices': list(p.indices),
                'normal': normal,
                'point': point,
                'vertices': [],
            }
        # loop all planes, find the intersection point of three plane
        keys = list(planes.keys())
        n = len(keys)
        vertices = []
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    point = threePlaneIntersection([planes[keys[i]],
                                                    planes[keys[j]],
                                                    planes[keys[k]]])
                    # remove point outside plane
                    point = convexhull(planes, point)
                    if point is not None:
                        planes[keys[i]]['vertices'].append(point)
                        planes[keys[j]]['vertices'].append(point)
                        planes[keys[k]]['vertices'].append(point)
        new_planes = {}
        for name, plane in planes.items():
            p = self.settings[plane['indices']]
            vertices, edges, faces = faces_from_vertices(
                plane['vertices'], plane['normal'])
            if len(vertices) >= 3:
                vertices += np.array(origin)
                new_planes[p.name] = self.get_plane_data(
                    vertices, edges, faces, p)
        self.crystal_planes = new_planes
        return new_planes

    def get_plane_data(self, vertices, edges, faces, p):
        """
        build edge
        """
        plane = {}
        if len(faces) > 0:
            plane = {'vertices': vertices,
                     'edges': edges,
                     'faces': faces,
                     'material_style': p.material_style,
                     'color': p.color,
                     'indices': p.indices,
                     'edges_cylinder': {'lengths': [], 'centers': [],
                                        'normals': [], 'vertices': 6,
                                        'color': (0.0, 0.0, 0.0, 1.0),
                                        'width': p.width,
                                        'battr_inputs': {},
                                        },
                     'battr_inputs': {'Bcrystalshape': p.as_dict()},
                     'show_edge': p.show_edge,
                     'slicing': p.slicing,
                     'boundary': p.boundary,
                     }
            for edge in edges:
                center = (vertices[edge[0]] + vertices[edge[1]])/2.0
                vec = vertices[edge[0]] - vertices[edge[1]]
                length = np.linalg.norm(vec)
                nvec = vec/length
                plane['edges_cylinder']['lengths'].append(length)
                plane['edges_cylinder']['centers'].append(center)
                plane['edges_cylinder']['normals'].append(nvec)
        return plane

    def draw(self, plane_name="ALL", no=None, origin=None):
        """Draw crystal shape
        no: int
            spacegroup of structure, if None, no will be determined by
            get_spacegroup_number()
        origin: xyz vector
            The center of cyrstal shape
        """
        from batoms.utils.butils import clean_coll_object_by_type
        # delete old plane
        clean_coll_object_by_type(self.batoms.coll, 'CRYSTALSHAPE')
        if no is not None:
            self.no = no
        if origin is None:
            origin = self.batoms.cell.origin
        planes = self.build_crystal(self.batoms.cell, origin=origin)
        for species, plane in planes.items():
            if plane_name.upper() != "ALL" and species != plane_name:
                continue
            name = '%s_%s_%s' % (self.label, self.name, species)
            self.delete_obj(name)
            obj = draw_surface_from_vertices(name, plane,
                                             coll=self.batoms.coll,
                                             )
            mat = self.build_materials(name, color=plane['color'],
                                       material_style=plane['material_style'])
            obj.data.materials.append(mat)
            obj.parent = self.batoms.obj
            obj.batoms.type = 'CRYSTALSHAPE'
            obj.batoms.label = self.batoms.label
            if plane['show_edge']:
                name = '%s_%s_%s' % (self.label, 'crystal_edge', species)
                self.delete_obj(name)
                obj = draw_cylinder(name=name,
                                    datas=plane['edges_cylinder'],
                                    coll=self.batoms.coll,
                                    )
                obj.parent = self.batoms.obj
                obj.batoms.type = 'CRYSTALSHAPE'
                obj.batoms.label = self.batoms.label

    @property
    def setting(self):
        from batoms.utils import deprecated
        """setting object."""
        deprecated('"setting" will be deprecated in the furture, please use "settings".')
        return self.settings

def save_image(data, filename, interpolation='bicubic'):
    """
    """
    import pylab as plt
    import numpy as np
    data = data.T
    size = np.shape(data)
    logger.debug('size: {}'.format(size))
    fig = plt.figure(figsize=(size[1]/size[0]*10, 10))
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(data, interpolation=interpolation)
    plt.savefig(filename, dpi=300)


def faces_from_vertices(vertices, normal,
                        include_center=False, scale=[1, 1, 1]):
    """
    get faces from vertices
    """
    # remove duplicative point
    vertices = np.unique(vertices, axis=0)
    n = len(vertices)
    if n < 3:
        return vertices, [], []
    center = np.mean(vertices, axis=0)
    v1 = vertices[0] - center
    angles = [[0, 0]]
    normal = normal/(np.linalg.norm(normal) + 1e-6)
    for i in range(1, n):
        v2 = vertices[i] - center
        x = np.cross(v1, v2)
        c = np.sign(np.dot(x, normal))
        angle = np.arctan2(c, np.dot(v1, v2))
        angles.append([i, angle])
    # scale
    vec = vertices - center
    # length = np.linalg.norm(vec, axis = 1)
    # nvec = vec/length[:, None]
    vertices = center + np.array([scale])*vec
    # search convex polyhedra
    angles = sorted(angles, key=lambda l: l[1])
    if not include_center:
        faces = [[a[0] for a in angles]]
        edges = [[angles[0][0], angles[-1][0]]]
        for i in range(0, n - 1):
            edges.append([angles[i][0], angles[i + 1][0]])
    else:
        # add center to vertices
        vertices = np.append(vertices, center.reshape(1, 3), axis=0)
        icenter = len(vertices) - 1
        faces = [[angles[0][0], angles[-1][0], icenter]]
        edges = [[angles[0][0], icenter],
                 [angles[0][0], angles[-1][0]],
                 [angles[-1][0], icenter]]
        for i in range(0, n - 1):
            faces.append([angles[i][0], angles[i + 1][0], icenter])
            edges.extend([[angles[i][0], angles[i + 1][0]],
                          [angles[i][0], icenter],
                          [angles[i + 1][0], icenter]])
    return vertices, edges, faces


def linePlaneIntersection(line, normal, point):
    """
    3D Line Segment and Plane Intersection
    - Point
    - Line contained in plane
    - No intersection
    """
    d = np.dot(point, normal)
    normalLine = line[0] - line[1]
    a = np.dot(normalLine, normal)
    # No intersection or Line contained in plane
    if np.isclose(a, 0):
        return None
    # in same side
    b = np.dot(line, normal) - d
    if b[0]*b[1] > 0:
        return None
    # Point
    v = point - line[0]
    d = np.dot(v, normal)/a
    point = np.round(line[0] + normalLine*d, 6)
    return point


def threePlaneIntersection(planes):
    """
    3D three Planes Intersection
    """
    # Condition for a single point of intersection
    mat = np.concatenate([[plane['normal']] for plane in planes], axis=0)
    det = np.linalg.det(mat)
    if np.isclose(det, 0):
        return None
    darray = np.array([np.dot(plane['point'], plane['normal'])
                      for plane in planes])
    point = np.zeros(3)
    for i in range(3):
        temp = mat.copy()
        temp[:, i] = darray
        point[i] = round(np.linalg.det(temp)/det, 6)
    return point


def convexhull(planes, point):
    """
    find points at the same side of origin for all planes
    Todo: '+' or '-'
    """
    if point is None:
        return None
    for name, plane in planes.items():
        x1 = np.dot(point, plane['normal']) - \
            np.dot(plane['point'], plane['normal'])
        x2 = np.dot(np.array([0, 0, 0]), plane['normal']) - \
            np.dot(plane['point'], plane['normal'])
        if abs(x1) > 1e-6 and x1*x2 < -1e-6:
            return None
    return point

    def as_dict(self):
        """
        """
        data = {}
        data['settings'] = self.settings.as_dict()
        data.update(self.settings.bpy_data.as_dict())
        return data
