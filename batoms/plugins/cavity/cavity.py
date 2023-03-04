"""Definition of the Cavity class.

This module defines the Cavity object in the Batoms package.

"""

import bpy
from time import time
import numpy as np
from batoms.attribute import Attributes
from batoms.base.object import ObjectGN
from batoms.plugins.base import PluginObject
from .setting import CavitySettings
from scipy import spatial
from batoms.utils.butils import object_mode, get_nodes_by_name
from batoms.utils import string2Number
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)

default_attributes = [
    {"name": 'species_index', "data_type": 'INT', "dimension": 0},
    {"name": 'species', "data_type": 'STRING', "dimension": 0},
    {"name": 'show', "data_type": 'INT', "dimension": 0},
    {"name": 'scale', "data_type": 'FLOAT', "dimension": 0},
]

default_GroupInput = [
    ['species_index', 'NodeSocketInt'],
    ['show', 'NodeSocketBool'],
    ['scale', 'NodeSocketFloat'],
]

default_cavity_datas = {
    'species_index': np.ones(0, dtype=int),
    'centers': np.zeros((1, 0, 3)),
    'show': np.zeros(0, dtype=int),
    'scale': np.ones(0, dtype=float),
}


class Cavity(ObjectGN, PluginObject):
    def __init__(self,
                 label=None,
                 cavity_datas=None,
                 location=np.array([0, 0, 0]),
                 batoms=None,
                 resolution = None,
                 minCave = None,
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
        ObjectGN.__init__(self, label, name)
        PluginObject.__init__(self, label, 'cavity', batoms)
        if resolution is not None:
            self.resolution = resolution
        if minCave is not None:
            self.minCave = minCave
        flag = self.load()
        if not flag:
            if cavity_datas is None:
                cavity_datas = default_cavity_datas
            self.build_object(cavity_datas)
            self.settings = CavitySettings(
                self.label, batoms=batoms, parent=self)
        else:
            self.settings = CavitySettings(
                self.label, batoms=batoms, parent=self)
            self._attributes = Attributes(label=self.label, parent=self, obj_name=self.obj_name)
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

    def draw(self):
        """Draw cavity.
        """
        # delete old cavity
        cavity_datas = self.build_cavity()
        self.set_arrays(cavity_datas)

    def build_kdtree(self, positions):
        """

        """
        from scipy import spatial
        tstart = time()
        self.kdtree = spatial.KDTree(positions)
        self.kdtree_mesh = spatial.KDTree(self.meshgrids)
        # print('KDTree positions: %s' % (time() - tstart))

    def query_distance(self, points, parallel=1):
        tstart = time()
        distance, indices = self.kdtree.query(points, workers=parallel)
        # print('KDTree query: %s' % (time() - tstart))
        # here we consider the size of atoms using one value for all.
        distance = distance - self.atomRadius
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

    def base_meshgrids(self, grids):
        x, y, z = np.meshgrid(grids[0], grids[1], grids[2],
                              indexing='ij')  # , sparse=True)
        meshgrids = np.c_[x.ravel(), y.ravel(), z.ravel()]
        return meshgrids


    def build_grid(self, cell, resolution):
        """
        generate gridpoints
        """
        length = cell.length
        npoint = [int(l/resolution) for l in length]
        grids = []
        for i in range(3):
            grids.append(np.arange(npoint[i])/npoint[i])
        meshgrids = self.base_meshgrids(grids)
        meshgrids = np.dot(meshgrids, cell.array)
        self.meshgrids = meshgrids
        self.shape = npoint

    def build_cavity(self):
        """Find cage base on radius range and boundary range
        Algorith:
            1) build a meshgrid
            2) calulated first neighbour distance by kdtree
            3) find the max distance (max radius sphere)
            4) remove meshgrid within this sphere,
            5) find next max radius, ..., untill radius < minCave
            6) remove all spheres contact with boundary
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
        spheres = self.find_cage_spheres(distances)
        spheres = self.check_sphere_boundary(spheres, cell)
        # spheres = self.refine_spheres(spheres)
        ns = len(spheres['centers'])
        color_names = list(basic_colors.keys())
        ic = 0
        # the default setting
        ncav = len(self.settings.bpy_setting)
        #
        indices = []
        for i in range(ns):
            has_r = False
            for cav in self.settings.bpy_setting:
                if spheres['radii'][i] > cav.min and spheres['radii'][i] < cav.max:
                    has_r = True
                    indices.append(i)
            # init setting if not exist (ncav == 0)
            if not has_r and ncav == 0:
                name = len(self.settings.bpy_setting)
                cav = {'min': np.floor(spheres['radii'][i]), 'max': np.ceil(
                    spheres['radii'][i]), 'color': basic_colors[color_names[ic]]}
                self.settings[name] = cav
                indices.append(i)
                ic += 1
                if ic == len(color_names):
                    ic = 0
        #
        cavities = {}
        for key, value in spheres.items():
            if len(spheres[key].shape)==2:
                cavities[key] = spheres[key][indices, :]
            else:
                cavities[key] = spheres[key][indices]
        nc = len(cavities['centers'])
        species_index = np.zeros(nc, dtype = int)
        show = np.ones(nc, dtype=int)
        scale = np.ones(nc, dtype=int)
        ncollection = len(self.settings.bpy_setting)
        for i in range(ncollection):
            cav = self.settings.bpy_setting[i]
            indices = np.where((cavities['radii'] > cav.min) & (
                cavities['radii'] < cav.max))[0]
            species_index[indices] = string2Number(cav.name)
            scale[indices] = cavities['radii'][indices]*cav.scale
        # logger.debug(scale)
        cavities.update({'centers': np.array([cavities['centers']]),
                        'species_index': species_index,
                        'scale': scale,
                        'show': show})
        return cavities

    def find_cage_spheres(self, distances):
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
        # refine center
        center, dmax = self.refine_spheres(center)
        centers = []
        radii = []
        while dmax > self.minCave:
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
            center, dmax = self.refine_spheres(center)
            # print(center, dmax)
        spheres = {'centers': np.array(centers), 'radii': np.array(radii)}
        # print('spheres: ', spheres)
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
        # print('spheres after check boundary: ', spheres)
        return spheres

    def refine_spheres(self, center):
        """Refine centers of spheres

        Find center off meshgrids.

        Args:
            spheres (_type_): _description_
        """
        # make dense meshgrid ground centers of spheres
        npoint = [5, 5, 5]
        grids = []
        for i in range(3):
            grid = np.arange(-npoint[i], npoint[i] + 1) / \
                npoint[i]*self.resolution/2
            grids.append(grid)
        meshgrids0 = self.base_meshgrids(grids)
        meshgrids = meshgrids0 + center
        indices, distances = self.query_distance(meshgrids)
        imax = np.argmax(distances)
        radius = distances[imax]
        return center, radius

    def build_object(self, cavity_datas, attributes = {}):
        """Build the main Batoms object

        Args:
            label (str):
                Name of the object
            arrays (array):
                arrays of properties for each atoms
            location (list, optional):
                Location of the object. Defaults to [0, 0, 0].
        """
        if len(cavity_datas['centers'].shape) == 2:
            self._frames = {'centers': np.array([cavity_datas['centers']]),
                            }
            centers = cavity_datas['centers']
        elif len(cavity_datas['centers'].shape) == 3:
            self._frames = {'centers': cavity_datas['centers'],
                            }
            centers = cavity_datas['centers'][0]
        else:
            raise Exception('Shape of centers is wrong!')
        attributes.update({
            'species_index': cavity_datas['species_index'],
            'show': cavity_datas['show'],
            'scale': cavity_datas['scale'],
            })
        name = self.obj_name
        self.delete_obj(name)
        mesh = bpy.data.meshes.new(name)
        obj = bpy.data.objects.new(name, mesh)
        obj.data.from_pydata(centers, [], [])
        obj.batoms.type = 'CAVITY'
        obj.batoms.label = self.batoms.label
        self.batoms.coll.objects.link(obj)
        self._attributes = Attributes(label=self.label, parent=self, obj_name=self.obj_name)
        # Add attributes
        self._attributes.add(default_attributes)
        # add cell object as its child
        obj.parent = self.batoms.obj
        self.set_attributes(attributes)
        self.build_geometry_node()

    def build_geometry_node(self):
        """Geometry node for instancing sphere on vertices!
        """
        from batoms.utils.butils import build_modifier
        name = 'GeometryNodes_%s' % self.obj_name
        modifier = build_modifier(self.obj, name)
        inputs = modifier.node_group.inputs
        GroupInput = modifier.node_group.nodes[0]
        GroupOutput = modifier.node_group.nodes[1]
        # add new output sockets
        for att in default_GroupInput:
            inputs.new(att[1], att[0])
            id = inputs[att[0]].identifier
            modifier['%s_use_attribute' % id] = True
            modifier['%s_attribute_name' % id] = att[0]
        gn = modifier
        # print(gn.name)
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_JoinGeometry' % self.label,
                                         'GeometryNodeJoinGeometry')
        SeparateGeometry = \
            get_nodes_by_name(gn.node_group.nodes,
                              '%s_SeparateGeometry' % self.label,
                              'GeometryNodeSeparateGeometry')
        gn.node_group.links.new(GroupInput.outputs['Geometry'],
                                SeparateGeometry.inputs['Geometry'])
        gn.node_group.links.new(GroupInput.outputs[2],
                                SeparateGeometry.inputs['Selection'])
        gn.node_group.links.new(SeparateGeometry.outputs[0],
                                JoinGeometry.inputs['Geometry'])
        gn.node_group.links.new(JoinGeometry.outputs['Geometry'],
                                GroupOutput.inputs['Geometry'])

    def add_geometry_node(self, spname, instancer):
        """Add geometry node for each species

        Args:
            spname (str):
                Name of the species
            instancer (bpy.data.object):
                Object to be instanced
        """
        from batoms.utils.butils import compareNodeType
        gn = self.gnodes
        GroupInput = gn.node_group.nodes[0]
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                                         '%s_JoinGeometry' % self.label,
                                         'GeometryNodeJoinGeometry')
        CompareSpecies = get_nodes_by_name(gn.node_group.nodes,
                                           'CompareFloats_%s_%s' % (
                                               self.label, spname),
                                           compareNodeType)
        CompareSpecies.operation = 'EQUAL'
        # CompareSpecies.data_type = 'INT'
        CompareSpecies.inputs[1].default_value = string2Number(spname)
        InstanceOnPoint = get_nodes_by_name(gn.node_group.nodes,
                                            'InstanceOnPoint_%s_%s' % (
                                                self.label, spname),
                                            'GeometryNodeInstanceOnPoints')
        ObjectInfo = get_nodes_by_name(gn.node_group.nodes,
                                       'ObjectInfo_%s_%s' % (
                                           self.label, spname),
                                       'GeometryNodeObjectInfo')
        ObjectInfo.inputs['Object'].default_value = instancer
        BoolShow = get_nodes_by_name(gn.node_group.nodes,
                                     'BooleanMath_%s_%s_1' % (
                                         self.label, spname),
                                     'FunctionNodeBooleanMath')
        #
        gn.node_group.links.new(GroupInput.outputs['Geometry'],
                                InstanceOnPoint.inputs['Points'])
        gn.node_group.links.new(GroupInput.outputs[1],
                                CompareSpecies.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[2], BoolShow.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[3],
                                InstanceOnPoint.inputs['Scale'])
        gn.node_group.links.new(CompareSpecies.outputs[0], BoolShow.inputs[1])
        gn.node_group.links.new(BoolShow.outputs['Boolean'],
                                InstanceOnPoint.inputs['Selection'])
        gn.node_group.links.new(ObjectInfo.outputs['Geometry'],
                                InstanceOnPoint.inputs['Instance'])
        gn.node_group.links.new(InstanceOnPoint.outputs['Instances'],
                                JoinGeometry.inputs['Geometry'])

    def update_geometry_node_instancer(self):
        """
        """
        for sp in self.settings.bpy_setting:
            self.settings.build_instancer(sp.as_dict())

    def set_arrays(self, arrays):
        """
        """
        # if len(arrays['positions']) == 0:
        #     return
        attributes = self.attributes
        # same length
        dnvert = len(arrays['species_index']) - \
            len(attributes['species_index'])
        if dnvert > 0:
            self.add_vertices_bmesh(dnvert)
        elif dnvert < 0:
            self.delete_vertices_bmesh(range(-dnvert))
        self.set_frames(arrays)
        self.set_attributes({'species_index': arrays['species_index']})
        self.set_attributes({'scale': arrays['scale']})
        self.set_attributes({'show': arrays['show']})
        self.update_mesh()
        self.update_geometry_node_instancer()

    def set_frames(self, frames=None, frame_start=0, only_basis=False):
        if frames is None:
            frames = self._frames
        nframe = len(frames['centers'])
        if nframe == 0:
            return
        name = '%s_cavity' % (self.label)
        obj = self.obj
        self.set_obj_frames(name, obj, frames['centers'])

    @property
    def objs(self):
        objs = {}
        name = '{}_{}'.format(self.label, self.name)
        obj = bpy.data.objects.get(name)
        objs[name] = obj
        return objs

    @property
    def setting(self):
        from batoms.utils import deprecated
        """setting object."""
        deprecated('"setting" will be deprecated in the furture, please use "settings".')
        return self.settings


    def as_dict(self):
        """
        """
        data = {}
        data['settings'] = self.settings.as_dict()
        data.update(self.settings.bpy_data.as_dict())
        return data

    @property
    def minCave(self):
        return self.batoms.coll.Bcavity.minCave

    @minCave.setter
    def minCave(self, minCave):
        self.batoms.coll.Bcavity.minCave = minCave

    @property
    def resolution(self):
        return self.batoms.coll.Bcavity.resolution

    @resolution.setter
    def resolution(self, resolution):
        self.batoms.coll.Bcavity.resolution = resolution

    @property
    def atomRadius(self):
        return self.batoms.coll.Bcavity.atomRadius

    @atomRadius.setter
    def atomRadius(self, atomRadius):
        self.batoms.coll.Bcavity.atomRadius = atomRadius
