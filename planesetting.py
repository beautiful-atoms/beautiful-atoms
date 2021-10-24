"""
"""
import bpy
from batoms.bondsetting import Setting
import numpy as np
from time import time

class PlaneSetting(Setting):
    """
    PlaneSetting object

    The PlaneSetting object store the polyhedra information.

    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.

    """
    def __init__(self, label, plane = None) -> None:
        Setting.__init__(self, label)
        self.name = 'bplane'
        if plane is not None:
            for key, data in plane.items():
                self[key] = data
    def add(self, indices):
        self[indices] = {'indices': indices}
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Indices   distance          color \n'
        for p in self.collection:
            s += '{0:10s}   {1:1.3f}   [{2:1.1f}  {3:1.1f}  {4:1.1f} {5:1.1f} ] \n'.format(\
                p.name, p.distance, p.color[0], p.color[1], p.color[2], p.color[3])
        s += '-'*60 + '\n'
        return s
    def build_crystal(self, no, bcell):
        """
        Build crystal

        no: Int
            spacegroup number
        """
        from scipy.spatial import ConvexHull
        from batoms.tools import get_equivalent_indices
        cellArray = bcell.array
        planes = []
        plane_kinds = {}
        for p in self:
            indices = get_equivalent_indices(no, p.indices)
            for i in indices:
                p1 = p.as_dict()
                normal = np.dot(i, cellArray)
                point = p.distance*normal
                p1.update({'name': '%s-%s-%s'%(i[0], i[1], i[2]), 
                            'indices': i, 
                            'normal': normal, 
                            'point':point, 
                            'vertices': [],
                            'edges': [],
                            'faces': [],
                            })
                planes.append(p1)
        # loop all planes, find the intersection point of three plane
        n = len(planes)
        vertices = []
        points = []
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    point = threePlaneIntersection([planes[i], planes[j], planes[k]])
                    if point is not None:
                        planes[i]['vertices'].append(point)
                        planes[j]['vertices'].append(point)
                        planes[k]['vertices'].append(point)
        plane_kinds = {}
        for plane in planes:
            vertices, edges, faces = faces_from_vertices(plane['vertices'], plane['normal'])
            if len(faces) > 0:
                plane_kinds[plane['name']] = {'vertices': vertices, 
                              'edges': edges, 
                              'faces': faces,
                              'color': plane['color'],
                              'indices': plane['indices'],
                              }
        return plane_kinds

    def build_plane(self, bcell):
        """
        Build vertices, edges and faces of plane.
        
        """
        from scipy.spatial import ConvexHull
        from batoms.tools import get_polyhedra_kind
        from scipy.spatial import ConvexHull
        planes = {}
        cellEdges = bcell.edges
        cellVerts = bcell.verts
        cellArray = bcell.array
        planes = {}
        for p in self:
            intersect_points = []
            normal = np.dot(p.indices, cellArray)
            # get intersection point
            for edge in cellEdges:
                line = cellVerts[edge]
                point = p.distance*normal
                intersect_point = linePlaneIntersection(line, normal, point)
                if intersect_point is not None:
                    intersect_points.append(intersect_point)
                # get verts, edges, faces by Hull
            vertices, edges, faces = faces_from_vertices(intersect_points, normal)
            if len(faces) > 0:
                planes[p.name] = {'vertices': vertices, 
                              'edges': edges, 
                              'faces': faces,
                              'color': p.color,
                              'indices': p.indices,
                              }
        self.planes = planes
        return planes

def faces_from_vertices(vertices, normal):
    """
    get faces from vertices
    """
    n = len(vertices)
    if n < 3: return vertices, [], []
    center = np.mean(vertices)
    v1 = vertices[0] - center
    angles = [[0, 0]]
    for i in range(1, n):
        v2 = vertices[i] - center
        x = np.cross(v1, v2)
        normal = normal/(np.linalg.norm(normal) + 1e-6)
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
    point = line[0] + normalLine*d
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
        point[i] = np.linalg.det(temp)/det
    return point
    