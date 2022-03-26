"""Definition of the cavity class.

This module defines the cavity object in the Batoms package.

"""

import bpy
from time import time
import numpy as np
from batoms.base.object import BaseObject
from batoms.cavity.cavitysetting import cavitySettings


class Cavity(BaseObject):
    def __init__(self,
                 label=None,
                 location=np.array([0, 0, 0]),
                 batoms=None,
                 ):
        """Cavity Class

        Args:
            label (_type_, optional): _description_. Defaults to None.
            location (_type_, optional): _description_. Defaults to np.array([0, 0, 0]).
            batoms (_type_, optional): _description_. Defaults to None.
        """
        #
        self.batoms = batoms
        self.label = label
        name = 'cavity'
        BaseObject.__init__(self, label, name)
        self.setting = cavitySettings(
            self.label, batoms=batoms, parent=self)

    def build_cavity(self, cell):
        volume = self.batoms.volume
        cavity = {}
        self.build_kdtree()
        self.build_grid()
        self.build_distances()
        for cav in self.setting.collection:
            name = cav.name
            min = cav.min
            max = cav.max
            color = cav.color
            verts, faces = self.find_cage_sphere(volume, cell, level)
            cavity[name] = {'vertices': verts,
                            'edges': [],
                            'faces': faces,
                            'color': color,
                            'battr_inputs': {'cavity': cav.as_dict()}
                            }
        return cavity

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

    def draw(self, cavity_name='ALL'):
        """Draw cavity.
        """
        from batoms.draw import draw_surface_from_vertices
        from batoms.utils.butils import clean_coll_object_by_type
        # delete old cavity
        clean_coll_object_by_type(self.batoms.coll, 'cavity')
        cavity = self.build_cavity(self.batoms.cell)
        for name, cavity_data in cavity.items():
            if cavity_name.upper() != "ALL" and name.name != cavity_name:
                continue
            name = '%s_%s_%s' % (self.label, 'cavity', name)
            self.delete_obj(name)
            obj = draw_surface_from_vertices(name,
                                             datas=cavity_data,
                                             coll=self.batoms.coll,
                                             )
            obj.batoms.type = 'cavity'
            obj.batoms.label = self.label
            obj.parent = self.batoms.obj
            # material
            mat = self.build_materials(name, cavity_data['color'])
            obj.data.materials.append(mat)

    def draw_cavity_sphere(self, radius, boundary=[[0, 1], [0, 1], [0, 1]]):
        """
        cavity
        for porous materials
        >>> from ase.io import read
        >>> atoms = read('docs/source/_static/datas/mof-5.cif')
        """
        from batoms.tools import find_cage_sphere
        object_mode()
        self.clean_atoms_objects('ghost')
        positions = find_cage_sphere(
            self.cell, self.atoms.positions, radius, boundary=boundary)
        ba = Batom(self.label, 'X', positions, scale=radius*0.9,
                   material_style='default', bsdf_inputs=self.bsdf_inputs,
                   color_style=self.color_style)

        # ba.color = [ba.color[0], ba.color[1], ba.color[2], 0.8]
        self.coll.children['%s_ghost' % self.label].objects.link(ba.batom)
        self.coll.children['%s_ghost' % self.label].objects.link(ba.instancer)

    def build_kdtree(self):
        """
        Algorithm:
        Use KDTree to find the nearest atom for all vertices with
        the power distance as metric.
        """
        from scipy import spatial
        tstart = time()
        self.kdtree = spatial.KDTree(self.batoms.positions)
        print('KDTree positions: %s' % (time() - tstart))
    
    def query_distance(self, points, parallel = 1):
        tstart = time()
        distance, indices = self.kdtree.query(points, workers=parallel)
        print('KDTree query: %s' % (time() - tstart))
        return indices, distance
    
    def query_radius(self, points, radii, k=1):
        """
        Algorithm:
        Use KDTree to query the tree for neighbors within a radius r.
        """
        from scipy import spatial
        tstart = time()
        # ----------------------------------------------------
        tstart = time()
        indices = self.kdtree.query_ball_point(points, radii)
        print('KDTree query: %s' % (time() - tstart))
        tstart = time()
        indices = np.unique(np.concatenate(indices).astype(int))
        print('array indices: %s' % (time() - tstart))
        return indices

    def build_grid(self, resolution):
        """
        generate gridpoints
        """
        from batoms.utils import build_grid
        cell = self.batoms.cell.arrays
        length = self.batoms.cell.length
        npoint = [int(l/resolution) for l in length]
        self.meshgrids, self.shape = build_grid([[0, 1], [0, 1], [0, 1]], resolution)
        grids = []
        for i in range(3):
            grids.append(np.arange(0, 1, npoint[i]))
        x, y, z = np.meshgrid(grids[0], grids[1], grids[2],
                            indexing='ij')  # , sparse=True)
        meshgrids = np.c_[x.ravel(), y.ravel(), z.ravel()]
        meshgrids = np.dot(meshgrids, cell)
        dis, indices = self.query(meshgrids)
        self.dis = dis
        self.indices = indices


    def find_cage_sphere(self, min, max):
        """Find cage base on radius range and boundary range

        Args:
            cell (_type_): _description_
            positions (_type_): _description_
            radius (_type_): _description_
            step (float, optional): _description_. Defaults to 1.0.
            boundary (list, optional): _description_. Defaults to [[0, 1], [0, 1], [0, 1]].

        Returns:
            _type_: _description_
        """
        # dists<3.0
        dis = self.dis
        imax = np.argmax(dis)
        dis[imax]
        indices = [[imax, dmax]]
        mask = []
        while dmax > min:
            mask1 = self.query_radius()
            mask.extend(mask1)
            dis = dis[~mask]
            imax = np.argmax(dis[~mask])
            dmax = dis[imax]
            indices.extend([imax, dmax])
        return indices