"""

Lattice Planes

To insert lattice planes in structural models.
"""
from numpy.core.records import array
import bpy
from batoms.bondsetting import Setting, tuple2string
import numpy as np
from time import time
from batoms.tools import get_equivalent_indices

class PlaneSetting(Setting):
    """
    PlaneSetting object

    The PlaneSetting object store the polyhedra information.

    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.

    """
    def __init__(self, label, no = None, plane = None) -> None:
        Setting.__init__(self, label)
        self.name = 'bplane'
        self.no = no
        if plane is not None:
            for key, data in plane.items():
                self[key] = data
    def __setitem__(self, index, setdict):
        """
        Set properties
        """
        name = tuple2string(index)
        p = self.find(name)
        if p is None:
            p = self.collection.add()
        p.indices = index
        p.name = name
        p.flag = True
        for key, value in setdict.items():
            setattr(p, key, value)
        p.label = self.label
        if p.symmetry:
            setdict = p.as_dict()
            indices = get_equivalent_indices(self.no, p.indices)
            for index in indices:
                name = tuple2string(index)
                p1 = self.find(name)
                if p1 is None:
                    p1 = self.collection.add()
                for key, value in setdict.items():
                    setattr(p1, key, value)
                p1.name = name
                p1.indices = index
                p1.label = self.label

    def add(self, indices):
        self[indices] = {'indices': indices}
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Indices   distance  crystal   symmetry \n'
        for p in self.collection:
            s += '{0:10s}   {1:1.3f}   {2:10s}  {3:10s} \n'.format(\
                p.name, p.distance, str(p.crystal), str(p.symmetry))
        s += '-'*60 + '\n'
        return s
    def get_symmetry_indices(self):
        if self.no is None: return
        for p in self:
            if p.symmetry:
                indices = get_equivalent_indices(self.no, p.indices)
                for index in indices:
                    name = tuple2string(index)
                    p1 = self.find(name)
                    if p1 is None:
                        p1 = self.collection.add()
                        setdict = p.as_dict()
                        for key, value in setdict.items():
                            setattr(p1, key, value)
                        p1.name = name
                        p1.indices = index
                        p1.label = self.label
                        p1.flag = True
    def build_crystal(self, bcell):
        """
        Build crystal

        no: Int
            spacegroup number
        """
        from scipy.spatial import ConvexHull
        self.get_symmetry_indices()
        cellArray = bcell.array
        planes = {}
        for p in self:
            if not p.crystal: continue
            normal = np.dot(p.indices, cellArray)
            normal = normal/np.linalg.norm(normal)
            point = p.distance*normal
            planes[p.name] = {
                        'indices': p.indices, 
                        'normal': normal, 
                        'point':point, 
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
            p = self[plane['indices']]
            vertices, edges, faces = faces_from_vertices(plane['vertices'], plane['normal'])
            if len(vertices) > 3:
                new_planes[p.name] = self.get_plane_data(vertices, edges, faces, p)
        return new_planes

    def build_plane(self, bcell):
        """
        Build vertices, edges and faces of plane.
        
        """
        from scipy.spatial import ConvexHull
        from batoms.tools import get_polyhedra_kind
        from scipy.spatial import ConvexHull
        self.get_symmetry_indices()
        cellEdges = bcell.edges
        cellVerts = bcell.verts
        cellArray = bcell.array
        planes = {}
        for p in self:
            if p.crystal: continue
            intersect_points = []
            normal = np.dot(p.indices, cellArray)
            normal = normal/np.linalg.norm(normal)
            # get intersection point
            for edge in cellEdges:
                line = cellVerts[edge]
                point = p.distance*normal
                intersect_point = linePlaneIntersection(line, normal, point)
                if intersect_point is not None:
                    intersect_points.append(intersect_point)
                # get verts, edges, faces by Hull
            if len(intersect_points) < 3: continue
            vertices, edges, faces = faces_from_vertices(intersect_points, normal)
            planes[p.name] = self.get_plane_data(vertices, edges, faces, p)
        self.planes = planes
        return planes
    def get_plane_data(self, vertices, edges, faces, p):
        """
        build edge
        """
        plane = {}
        if len(faces) > 0:
            plane= {'vertices': vertices, 
                    'faces': faces,
                    'color': p.color,
                    'indices': p.indices,
                    'edges': {'lengths': [], 'centers': [], 
                            'normals': [], 'vertices': 6,
                            'color': (0.0, 0.0, 0.0, 1.0),
                            'width': p.width,
                            'battr_inputs': {},
                            }, 
                    'battr_inputs': {'bplane': p.as_dict()},
                    'show_edge': p.show_edge,
                    }
            for edge in edges:
                center = (vertices[edge[0]] + vertices[edge[1]])/2.0
                vec = vertices[edge[0]] - vertices[edge[1]]
                length = np.linalg.norm(vec)
                nvec = vec/length
                plane['edges']['lengths'].append(length)
                plane['edges']['centers'].append(center)
                plane['edges']['normals'].append(nvec)
        return plane
def faces_from_vertices(vertices, normal):
    """
    get faces from vertices
    """
    from scipy.spatial.distance import cdist
    # remove duplicative point
    vertices = np.unique(vertices, axis = 0)
    n = len(vertices)
    if n < 3: return vertices, [], []
    center = np.mean(vertices, axis = 0)
    v1 = vertices[0] - center
    angles = [[0, 0]]
    normal = normal/(np.linalg.norm(normal) + 1e-6)
    for i in range(1, n):
        v2 = vertices[i] - center
        x = np.cross(v1, v2)
        c = np.sign(np.dot(x, normal))
        angle = np.arctan2(c, np.dot(v1, v2))
        angles.append([i, angle])
    # search convex polyhedra
    angles = sorted(angles,key=lambda l:l[1])
    faces = [[a[0] for a in angles]]
    edges = [[angles[0][0], angles[-1][0]]]
    for i in range(0, n - 1):
        edges.append([angles[i][0], angles[i + 1][0]])
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
    mat = np.concatenate([[plane['normal']] for plane in planes], axis = 0)
    det = np.linalg.det(mat)
    if np.isclose(det, 0):
        return None
    darray = np.array([np.dot(plane['point'], plane['normal']) for plane in planes])
    point = np.zeros(3)
    for i in range(3):
        temp = mat.copy()
        temp[:, i] = darray
        point[i] = round(np.linalg.det(temp)/det, 6)
    return point
def convexhull(planes, point):
    """
    find points at the same side of origin for all planes
    """
    Flag = True
    if point is None: return None
    for name, plane in planes.items():
        x1 = np.dot(point, plane['normal']) + np.dot(plane['point'], plane['normal'])
        x2 = np.dot(np.array([0, 0, 0]), plane['normal']) + np.dot(plane['point'], plane['normal'])
        if x1*x2 < -1e-6:
            flag = False
            return None
    return point

