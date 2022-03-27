"""Definition of the cavity class.

This module defines the cavity object in the Batoms package.

"""

import bpy
from time import time
import numpy as np
from batoms.base.object import BaseObject
from batoms.cavity.cavitysetting import CavitySettings
from scipy import spatial


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
        self.resolution = 0.5
        self.setting = CavitySettings(
            self.label, batoms=batoms, parent=self)

    
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
        cavities = self.build_cavity()
        for name, data in cavities.items():
            if cavity_name.upper() != "ALL" and name.name != cavity_name:
                continue
            name = '%s_%s_%s' % (self.label, 'cavity', name)
            mat = self.build_materials(name, data['color'])
            for i in range(len(data['centers'])):
                name1 = "%s_%s"%(name, i)
                self.delete_obj(name1)
                bpy.ops.mesh.primitive_uv_sphere_add(location = data['centers'][i],
                                                    radius=data['radii'][i])
                obj = bpy.context.view_layer.objects.active
                obj.name = name1
                obj.data.name = name1
                obj.batoms.atom.radius = data['radii'][i]
                obj.batoms.type = 'CAVITY'
                obj.batoms.label = self.label
                bpy.ops.object.shade_smooth()
                obj.parent = self.batoms.obj
                # material
                obj.data.materials.append(mat)

    def draw_cavity_sphere(self, radius, boundary=[[0, 1], [0, 1], [0, 1]]):
        """
        cavity
        for porous materials
        >>> from ase.io import read
        >>> atoms = read('docs/source/_static/datas/mof-5.cif')
        """
        from batoms.tools import find_cage_spheres_2
        object_mode()
        self.clean_atoms_objects('ghost')
        positions = find_cage_spheres_2(
            self.cell, self.atoms.positions, radius, boundary=boundary)
        ba = Batom(self.label, 'X', positions, scale=radius*0.9,
                   material_style='default', bsdf_inputs=self.bsdf_inputs,
                   color_style=self.color_style)

        # ba.color = [ba.color[0], ba.color[1], ba.color[2], 0.8]
        self.coll.children['%s_ghost' % self.label].objects.link(ba.batom)
        self.coll.children['%s_ghost' % self.label].objects.link(ba.instancer)

    def build_kdtree(self, positions):
        """
        
        """
        from scipy import spatial
        tstart = time()
        self.kdtree = spatial.KDTree(positions)
        self.kdtree_mesh = spatial.KDTree(self.meshgrids)
        print('KDTree positions: %s' % (time() - tstart))
    
    def query_distance(self, points, parallel = 1):
        tstart = time()
        distance, indices = self.kdtree.query(points, workers=parallel)
        print('KDTree query: %s' % (time() - tstart))
        return indices, distance
    
    def query_radius(self, meshgrids, points, radii, k=1):
        """
        Algorithm:
        Use KDTree to query the tree for neighbors within a radius r.
        """
        tstart = time()
        kdtree_mesh = spatial.KDTree(meshgrids)
        indices = kdtree_mesh.query_ball_point(points, radii)
        # print('KDTree query: %s' % (time() - tstart))
        return indices[0]

    def build_grid(self, cell, resolution):
        """
        generate gridpoints
        """
        length = cell.length
        npoint = [int(l/resolution) for l in length]
        grids = []
        for i in range(3):
            grids.append(np.arange(npoint[i])/npoint[i])
        x, y, z = np.meshgrid(grids[0], grids[1], grids[2],
                            indexing='ij')  # , sparse=True)
        meshgrids = np.c_[x.ravel(), y.ravel(), z.ravel()]
        meshgrids = np.dot(meshgrids, cell.array)
        self.meshgrids = meshgrids
        self.shape = npoint
    

    def build_cavity(self):
        """Find cage base on radius range and boundary range
        Algorith:
            1) build a meshgrid
            2) calulated first neighbour distance by kdtree
            3) find indices that d > min
            4) reshape to indices (3d)
            5) label by scipy.ndimage
            6) check components (support periodic)
            7) find minmum radius
        Args:
            cell (_type_): _description_
            positions (_type_): _description_
            radius (_type_): _description_
            step (float, optional): _description_. Defaults to 1.0.
            boundary (list, optional): _description_. Defaults to [[0, 1], [0, 1], [0, 1]].

        Returns:
            _type_: _description_
        """
        from batoms.data import basic_colors
        arrays = self.batoms.arrays
        cell = self.batoms.cell
        self.build_grid(cell, self.resolution)
        self.build_kdtree(arrays["positions"])
        indices, distances = self.query_distance(self.meshgrids)
        spheres = self.find_cage_spheres(distances, minRadius=5)
        spheres = self.check_sphere_boundary(spheres, cell)
        cavities = {}
        # init setting if not exist
        color_names = list(basic_colors.keys())
        ic = 0
        for r in spheres['radii']:
            has_r = False
            for cav in self.setting.collection:
                if r > cav.min and r < cav.max:
                    has_r = True
            if not has_r:
                cav = {'min': np.floor(r), 'max': np.ceil(r), 'color': basic_colors[color_names[ic]]}
                self.setting['%s_%s'%(cav['min'], cav['max'])] = cav
                ic += 1
                if ic == len(color_names):
                    ic =0
            
        for cav in self.setting.collection:
            indices = np.where((spheres['radii'] > cav.min) & (spheres['radii'] < cav.max))[0]
            cavities[cav.name] = {'centers': spheres['centers'][indices],
                                'radii': spheres['radii'][indices],
                                'color': cav.color,
                                # 'battr_inputs': {'cavities': cav.as_dict()}
                                'min': cav.min,
                                'max': cav.max,
                                }
        # print(cavities)
        return cavities

    def find_cage_spheres(self, distances, minRadius):
        """
        Loop all sphere > min

        Args:
            distances (_type_): _description_

        Returns:
            _type_: _description_
        """
        meshgrids = self.meshgrids
        imax = np.argmax(distances)
        dmax = distances[imax]
        center = meshgrids[imax]
        centers = []
        radii = []
        while dmax > minRadius:
            centers.append(center)
            radii.append(dmax)
            n = len(distances)
            mask = np.ones(n, dtype=bool)
            indices1 = self.query_radius(meshgrids, [center], [dmax])
            mask[imax] = False
            mask[indices1] = False
            distances = distances[mask]
            meshgrids = meshgrids[mask]
            imax = np.argmax(distances)
            dmax = distances[imax]
            center = meshgrids[imax]
            # print(center, dmax)
        spheres = {'centers': np.array(centers), 'radii': np.array(radii)}
        print('spheres: ', spheres)
        return spheres

    def check_sphere_boundary(self, spheres0, cell):
        """Remove sphere contact with boundary
        distance to cell < radius
        """
        from batoms.neighborlist import pointCellDistance
        centers = []
        radii = []
        if len(spheres0['centers']) == 0:
            spheres = {'centers': np.array(centers), 'radii': np.array(radii)}
            return spheres

        dis = pointCellDistance(spheres0['centers'], cell)
        # print('pointCellDistance: ', d)
        
        ns = len(spheres0['centers'])
        for i in range(ns):
            dmin = np.min(dis[:, :, i])
            # print(dmin, spheres0['radii'][i])
            if dmin > spheres0['radii'][i]:
                centers.append(spheres0['centers'][i])
                radii.append(spheres0['radii'][i])
        spheres = {'centers': np.array(centers), 'radii': np.array(radii)}
        print('spheres after check: ', spheres)
        return spheres


    def build_cavity_2(self):
        """Find cage base on radius range and boundary range
        Algorith:
            1) build a meshgrid
            2) calulated first neighbour distance by kdtree
            3) find indices that d > min
            4) reshape to indices (3d)
            5) label by scipy.ndimage
            6) check components (support periodic)
            7) find minmum radius
        Args:
            cell (_type_): _description_
            positions (_type_): _description_
            radius (_type_): _description_
            step (float, optional): _description_. Defaults to 1.0.
            boundary (list, optional): _description_. Defaults to [[0, 1], [0, 1], [0, 1]].

        Returns:
            _type_: _description_
        """
        arrays = self.batoms.arrays
        cell = self.batoms.cell
        self.build_kdtree(arrays["positions"])
        self.build_grid(cell, self.resolution)
        indices, distance = self.query_distance(self.meshgrids)
        cavities = {}
        for cav in self.setting.collection:
            cavities[cav.name] = {'centers': [],
                                'radii': [],
                                'color': cav.color,
                                # 'battr_inputs': {'cavities': cav.as_dict()}
                                'min': cav.min,
                                'max': cav.max,
                                }
            s = self.find_cage_spheres_2(distance, radius=(cav.min + cav.max)/2)
            cavities['centers'].extend(s['center'])
            cavities['radii'].extend(s['radius'])
        return cavities

    def find_cage_spheres_2(self, distance, radius = 5):
        """
         Algorith:
            1) build a meshgrid
            2) calulated first neighbour distance by kdtree
            3) find indices that d > min (e.g. 3)
            4) reshape to indices (3d)
            5) label by scipy.ndimage
            6) check components (support periodic)
            7) find minmum radius
        """
        from scipy import ndimage
        dis = distance.reshape(self.shape)
        meshgrids = self.meshgrids.reshape((self.shape[0], self.shape[1], self.shape[2], 3))
        mask = dis > radius
        # check connected poits
        labels, nb = ndimage.label(mask)
        print(nb)
        non_boundary_labels = self.find_cage_pbc_2(labels)
        spheres = []
        for l in non_boundary_labels:
            mask = labels == l
            points = meshgrids[mask, :]
            center = np.mean(points, axis = 0)
            d = spatial.distance.cdist([center], points)
            dmax = np.max(d) + radius
            print(l, center, dmax)
            spheres.extend({'label': l, 'center': center, 'radius': dmax})
        print(spheres)
        return spheres


    def find_cage_pbc_2(self, labels):
        """_summary_

        Args:
            labels (_type_): _description_
            nb (_type_): _description_
        """
        from scipy.sparse import csgraph, csr_matrix
        neighbourList = []
        tstart = time()
        for i in range(labels.shape[0]):
            for j in range(labels.shape[1]):
                if labels[i, j, 0] > 0 and labels[i, j, -1] > 0:
                    neighbourList.append([labels[i, j, 0], labels[i, j, -1], 0, 0, -1])
        for i in range(labels.shape[0]):
            for j in range(labels.shape[2]):
                if labels[i, 0, j] > 0 and labels[i, -1, j] > 0:
                    neighbourList.append([labels[i, 0, j], labels[i, -1, j], 0, -1, 0])
        for i in range(labels.shape[1]):
            for j in range(labels.shape[2]):
                if labels[0, i, j] > 0 and labels[-1, i, j] > 0:
                    neighbourList.append([labels[0, i, j], labels[-1, i, j], -1, 0, 0])
        neighbourList = np.array(neighbourList)
        neighbourList = np.unique(neighbourList, axis = 0)
        print(neighbourList)
        nlabel = np.max(labels) + 1
        non_boundary_labels = list(set(range(1, nlabel)).difference(set(neighbourList[:, 0]).union(set(neighbourList[:, 1]))))
        print("non_boundary_labels: ", non_boundary_labels)
        # find connected component
        # cageDatas = {}
        # ai = neighbourList[:, 0]
        # aj = neighbourList[:, 1]
        # nlist = len(neighbourList)
        # data = np.ones(nlist, dtype=int)
        # matrix = csr_matrix((data, (ai, aj)), shape=(nlabel, nlabel))
        # n_components, component_list = csgraph.connected_components(matrix)
        # print(n_components)
        # print(component_list)
        # for i in range(n_components):
        #     indices = np.where(component_list == i)[0]
        #     n = len(indices)
        #     if n < 2:
        #         continue
        #     cageDatas[i] = {'sub': []}
        #     cageDatas[i]['indices'] = indices
        #     cageDatas[i]['offsets'] = indices
        print('find_cage_pbc_2: %s' % (time() - tstart))
        return non_boundary_labels
        
