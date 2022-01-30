"""Definition of the polyhedras class.

This module defines the polyhedras object in the Batoms package.

"""

from pandas import offsets
from sympy import false, im
import bpy
import bmesh
from time import time
from batoms.butils import object_mode
from batoms.tools import string2Number
import numpy as np
from batoms.base import BaseObject
from batoms.polyhedrasetting import PolyhedraSettings

default_attributes = [
            ['atoms_index1', 'INT', 'POINT'],
            ['atoms_index2', 'INT', 'POINT'],
            ['species_index', 'INT', 'POINT'],
            ['face_species_index', 'INT', 'FACE'],
            ['show', 'BOOLEAN', 'POINT'],
            ['style', 'INT', 'POINT'],
        ]

default_polyhedra_datas = {
        'atoms_index1': np.ones(0, dtype = int),
        'atoms_index2': np.ones(0, dtype = int),
        'species_index': np.ones(0, dtype = int),
        'face_species_index': np.ones(0, dtype = int),
        'vertices':np.zeros((0, 3)),
        # 'vectors':np.zeros((0, 3)),
        'offsets':np.zeros((0, 3)),
        # 'eulers':np.eye(3),
        # 'lengths':np.zeros((0, 3)),
        'edges': [],
        'faces': [],
        'widths': np.ones(0, dtype = float),
        'show': np.zeros(0, dtype = int),
        'styles': np.zeros(0, dtype = int),
        # 'model_styles':np.ones(0, dtype = int),
        }

class Polyhedras(BaseObject):
    """Polyhedras Class
    
    A Polyhedras object is linked to this main collection in Blender. 

    Parameters:

    label: str
        Name of the Bpolyhedras.
    species: str
        species of the atoms.
    positions: array
        positions
    locations: array
        The object’s origin location in global coordinates.
    element: str or list
        element of the atoms, list for fractional Occupancy
    segments: list of 2 Int
        Number of segments used to draw the UV_Sphere
        Default: [32, 16]
    subdivisions: Int
        Number of subdivision used to draw the ICO_Sphere
        Default: 2
    color_style: str
        "JMOL", "ASE", "VESTA"
    radii_style: str
        "covelent", "vdw", "ionic"
    shape: Int
        0, 1, or 2. ["UV_SPHERE", "ICO_SPHERE", "CUBE"]

    Examples:

    >>> from batoms.bond import Bbond
    >>> c = Bbond('C', [[0, 0, 0], [1.2, 0, 0]])

    """
    

    
    def __init__(self, 
                label = None,
                polyhedra_datas = None,
                location = np.array([0, 0, 0]),
                batoms = None,
                 ):
        #
        self.batoms = batoms
        self.label = label
        obj_name = '%s_polyhedra_center'%(self.label)
        bobj_name = 'bbond'
        BaseObject.__init__(self, obj_name = obj_name, bobj_name = bobj_name)
        self.setting = PolyhedraSettings(self.label, batoms = batoms, polyhedras = self)
        flag = self.load()
        if not flag and polyhedra_datas is not None:
            self.build_object(polyhedra_datas)
    
    def build_object(self, polyhedra_datas, attributes = {}):
        object_mode()
        """
        build child object and add it to main objects.
        """
        tstart = time()
        if len(polyhedra_datas['vertices'].shape) == 2:
            self._frames = (np.array([polyhedra_datas['vertices']]), 
                            np.array([polyhedra_datas['offsets']]),
                            )
            vertices = polyhedra_datas['vertices']
            offsets = polyhedra_datas['offsets']
        elif len(polyhedra_datas['vertices'].shape) == 3:
            self._frames = (polyhedra_datas['vertices'], 
                            polyhedra_datas['offsets'], 
                            )
            vertices = polyhedra_datas['vertices'][0]
            offsets = polyhedra_datas['offsets'][0]
        else:
            raise Exception('Shape of vertices is wrong!')
        nbond = len(polyhedra_datas['vertices'])
        show = np.ones(nbond, dtype = int)
        attributes.update({
                            'atoms_index1': polyhedra_datas['atoms_index1'], 
                            'atoms_index2': polyhedra_datas['atoms_index2'], 
                            'species_index': polyhedra_datas['species_index'], 
                            'face_species_index': polyhedra_datas['face_species_index'], 
                            'show': show,
                            # 'model_style': polyhedra_datas['model_styles'],
                            'style': polyhedra_datas['styles'],
                            })
        name = '%s_polyhedra_center'%self.label
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink = True)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(vertices, polyhedra_datas['edges'], polyhedra_datas['faces'])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        for attribute in default_attributes:
            mesh.attributes.new(name = attribute[0], type = attribute[1], domain = attribute[2])
        self.setting.coll.objects.link(obj)
        #
        name = '%s_polyhedra_offset'%self.label
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.data.objects.remove(obj, do_unlink = True)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(offsets, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        self.setting.coll.objects.link(obj)
        obj.hide_set(True)
        bpy.context.view_layer.update()
        self.set_attributes(attributes)
        self.build_geometry_node()
        self.set_frames(self._frames, only_basis = True)
        self.assign_materials()
        print('polyhedras: build_object: {0:10.2f} s'.format(time() - tstart))
    
    def assign_materials(self):
        # sort element by occu
        me = self.obj.data
        me.materials.clear()
        for sp in self.setting:
            me.materials.append(self.setting.materials[sp.species][0])
       
    def load(self):
        flag = True
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            flag = False
        return flag
    
    @property
    def gnodes(self):
        return self.get_gnodes()
    
    @gnodes.setter
    def gnodes(self, gnodes):
        self.set_gnodes(gnodes)
    
    def get_gnodes(self):
        name = 'GeometryNodes_%s_polyhedra'%self.label
        modifier = self.obj.modifiers.get(name)
        if modifier is None:
            self.build_geometry_node()
        return modifier
    
    def set_gnodes(self, gnodes):
        pass

    def build_geometry_node(self):
        """
        Geometry node for everything!
        """
        from batoms.butils import get_nodes_by_name
        tstart = time()
        name = 'GeometryNodes_%s_polyhedra'%self.label
        modifier = self.obj.modifiers.new(name = name, type = 'NODES')
        modifier.node_group.name = name
        #------------------------------------------------------------------
        # select attributes
        GroupInput = modifier.node_group.nodes.get('Group Input')
        # add new output sockets
        for att in default_attributes:
            GroupInput.outputs.new(type = att[1], name = att[0])
        # the above codes not works. maybe bug in blender, 
        # we add this, maybe deleted in the future
        for i in range(1, 6):
            test = get_nodes_by_name(modifier.node_group.nodes, 
                            'BooleanMath_%s'%i,
                            'FunctionNodeCompareFloats')
            modifier.node_group.links.new(GroupInput.outputs[i], test.inputs[0])
        #
        i = 2
        for att in default_attributes:
            modifier['Input_%s_use_attribute'%i] = 1
            modifier['Input_%s_attribute_name'%i] = att[0]
            i += 1
        gn = modifier
        #------------------------------------------------------------------
        GroupOutput = gn.node_group.nodes.get('Group Output')
        #------------------------------------------------------------------
        # calculate bond vector, length, rotation based on the index
        # Get four positions from batoms, bond and the second bond for high order bond plane
        ObjectBatoms = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_ObjectBatoms'%self.label,
                    'GeometryNodeObjectInfo')
        ObjectBatoms.inputs['Object'].default_value = self.batoms.obj
        PositionBatoms = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_PositionBatoms'%(self.label),
                        'GeometryNodeInputPosition')
        TransferBatoms = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_TransferBatoms'%(self.label),
                    'GeometryNodeAttributeTransfer')
        TransferBatoms.mapping = 'INDEX'
        TransferBatoms.data_type = 'FLOAT_VECTOR'
        gn.node_group.links.new(ObjectBatoms.outputs['Geometry'], TransferBatoms.inputs['Target'])
        gn.node_group.links.new(PositionBatoms.outputs['Position'], TransferBatoms.inputs['Attribute'])
        gn.node_group.links.new(GroupInput.outputs[2], TransferBatoms.inputs['Index'])
        #------------------------------------------------------------------
        # add positions with offsets
        # transfer offsets from object self.obj_o
        ObjectOffsets = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_ObjectOffsets'%(self.label),
                        'GeometryNodeObjectInfo')
        ObjectOffsets.inputs['Object'].default_value = self.obj_o
        PositionOffsets = get_nodes_by_name(gn.node_group.nodes, 
                        '%s_PositionOffsets'%(self.label),
                        'GeometryNodeInputPosition')
        TransferOffsets = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_TransferOffsets'%self.label,
                    'GeometryNodeAttributeTransfer')
        TransferOffsets.mapping = 'INDEX'
        TransferOffsets.data_type = 'FLOAT_VECTOR'
        gn.node_group.links.new(ObjectOffsets.outputs['Geometry'], TransferOffsets.inputs['Target'])
        gn.node_group.links.new(PositionOffsets.outputs['Position'], TransferOffsets.inputs['Attribute'])
        # we need one add operation to get the positions with offset
        VectorAdd = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_VectorAdd'%(self.label),
                    'ShaderNodeVectorMath')
        VectorAdd.operation = 'ADD'
        gn.node_group.links.new(TransferBatoms.outputs[0], VectorAdd.inputs[0])
        gn.node_group.links.new(TransferOffsets.outputs[0], VectorAdd.inputs[1])
        # set positions
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                        '%s_SetPosition'%self.label, 
                        'GeometryNodeSetPosition')
        gn.node_group.links.new(GroupInput.outputs['Geometry'], SetPosition.inputs['Geometry'])
        gn.node_group.links.new(VectorAdd.outputs[0], SetPosition.inputs['Position'])
        gn.node_group.links.new(SetPosition.outputs['Geometry'], GroupOutput.inputs['Geometry'])
        
        # find kinds by the names of species
        for sp in self.setting:
            self.add_geometry_node(sp.as_dict())
        print('Build geometry nodes for polyhedras: %s'%(time() - tstart))

    def add_geometry_node(self, sp):
        from batoms.butils import get_nodes_by_name
        gn = self.gnodes
        GroupInput = gn.node_group.nodes.get('Group Input')
        GroupOutput = gn.node_group.nodes.get('Group Output')
        previousNode = GroupOutput.inputs['Geometry'].links[0].from_socket
        # print(previousNode)
        # we need two compares for one species,
            # vertices
        CompareSpecies = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_CompareSpecies_%s_vertex'%(self.label, sp["species"]),
                    'FunctionNodeCompareFloats')
        CompareSpecies.operation = 'EQUAL'
        CompareSpecies.inputs[1].default_value = string2Number(sp["species"])
        gn.node_group.links.new(GroupInput.outputs[3], CompareSpecies.inputs[0])
        # face
        CompareSpeciesFace = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_CompareSpecies_%s_face'%(self.label, sp["species"]),
                    'FunctionNodeCompareFloats')
        CompareSpeciesFace.operation = 'EQUAL'
        CompareSpeciesFace.inputs[1].default_value = string2Number(sp["species"])
        gn.node_group.links.new(GroupInput.outputs[4], CompareSpeciesFace.inputs[0])
        #
        setMaterialIndex = get_nodes_by_name(gn.node_group.nodes, 
                '%s_setMaterialIndex_%s'%(self.label, sp["species"]),
                    'GeometryNodeSetMaterialIndex')
        setMaterialIndex.inputs[2].default_value = self.setting.materials[sp["species"]][1]
        gn.node_group.links.new(previousNode, setMaterialIndex.inputs['Geometry'])
        gn.node_group.links.new(CompareSpeciesFace.outputs[0], setMaterialIndex.inputs['Selection'])        
        gn.node_group.links.new(setMaterialIndex.outputs['Geometry'], GroupOutput.inputs['Geometry'])

    def update_polyhedras(self, ):
        """Draw polyhedras.
        calculate bond in all farmes, and save
        get the max number of polyhedras for each pair
        draw the polyhedras
        add shape key
        the extra polyhedras use the find bond data.
        """
        
        from batoms.butils import clean_coll_objects
        # if not self.bondlist:
        object_mode()
        # clean_coll_objects(self.coll, 'bond')
        frames = self.batoms.get_frames()
        arrays = self.batoms.arrays
        size = arrays['radius']*arrays['scale']
        species = arrays['species']
        # frames_boundary = self.batoms.get_frames(self.batoms.batoms_boundary)
        # frames_search = self.batoms.get_frames(self.batoms.batoms_search)
        nframe = len(frames)
        polyhedra_datas = {}
        tstart = time()
        for f in range(nframe):
            print('update polyhedra: ', f)
            positions = frames[f]
            # if len(frames_boundary) > 0:
            #     positions_boundary = frames_boundary[f]
            #     positions = positions + positions_boundary
            # if len(frames_search) > 0:
            #     positions_search = frames_search[f]
            #     positions = positions + positions_search
            bondlist = self.bondlists
            if len(bondlist['positions']) == 0:
                self.batoms.bonds.update_bonds()
                bondlist = self.batoms.bonds.arrays
            polyhedra_kinds = self.calc_polyhedra_data(bondlist, positions, self.batoms.cell)
            if f == 0:
                polyhedra_datas = polyhedra_kinds
            
        # print('calc polyhedra: {0:10.2f} s'.format(time() - tstart))
            # nframe = len(polyhedra_datas['positions'])
            # if nframe == 0: continue
            # change to same length
            # find max
            # nbs = [polyhedra_data['positions'][i].shape[0] for i in range(nframe)]
            # nb_max = max(nbs)
            # frames_polyhedra = np.zeros((nframe, nb_max, 7))
            # for i in range(nframe):
            #     frames_polyhedra[i, 0:nbs[i], :] = polyhedra_data['positions'][i]
            #     frames_polyhedra[i, nbs[i]:, 0:3] = frames[i][0]
        if len(polyhedra_datas) == 0:
            return
        self.set_polyhedra_datas(polyhedra_datas)
        # self.coll.objects.link(bb.obj)
        # bpy.data.collections['Collection'].objects.unlink(bb.obj)
        # bb.set_frames()
        # bpy.context.scene.frame_set(self.batoms.nframe)
        print('draw polyhedra: {0:10.2f} s'.format(time() - tstart))

    @property
    def obj_o(self):
        return self.get_obj_o()
    
    def get_obj_o(self):
        name = '%s_polyhedra_offset'%self.label
        obj_o = bpy.data.objects.get(name)
        if obj_o is None:
            raise KeyError('%s object is not exist.'%name)
        return obj_o

    @property
    def arrays(self):
        return self.get_arrays()
    
    @arrays.setter
    def arrays(self, arrays):
        self.set_arrays(arrays)

    def get_arrays(self):
        """
        """
        object_mode()
        tstart = time()
        arrays = self.attributes
        arrays.update({'positions': self.positions[0],
                        'offsets': self.positions[1],
                        })
        # print('get_arrays: %s'%(time() - tstart))
        return arrays
    
    @property
    def bondlists(self):
        return self.get_bondlists()
    
    def get_bondlists(self):
        try:
            bondlists = self.batoms.bonds.arrays
        except:
            bondlists = None
        return bondlists
        
    @property
    def polyhedra_datas(self):
        return self.get_polyhedra_datas()
    
    @polyhedra_datas.setter    
    def polyhedra_datas(self, polyhedra_datas):
        self.set_polyhedra_datas(polyhedra_datas)
    
    def get_polyhedra_datas(self):
        return np.array(self.mesh.polyhedra_datas)
    
    def set_polyhedra_datas(self, polyhedra_datas):
        """
        """
        if len(polyhedra_datas['vertices']) == 0:
            return
        attributes = self.attributes
        # same length
        if len(polyhedra_datas['vertices']) == len(attributes['show']):
            self.set_positions(polyhedra_datas['vertices'],
                                polyhedra_datas['offsets'],
                                )
            self.set_attributes({'atoms_index1': polyhedra_datas['atoms_index1']})
            self.set_attributes({'atoms_index2': polyhedra_datas['atoms_index2']})
            self.set_attributes({'species_index': polyhedra_datas['species_index']})
            self.set_attributes({'face_species_index': polyhedra_datas['face_species_index']})
            # self.set_attributes({'order': polyhedra_datas['orders']})
            self.set_attributes({'style': polyhedra_datas['styles']})
        else:
            # add or remove vertices
            self.build_object(polyhedra_datas)

    @property
    def local_positions(self):
        return self.get_local_positions()
    
    def get_local_positions(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        n = len(self)
        local_positions = []
        for obj in [self.obj, self.obj_o]:
            co = np.empty(n*3, dtype=np.float64)
            obj.data.vertices.foreach_get('co', co)  
            local_positions.append(co.reshape((n, 3)))
        return local_positions
    
    @property
    def positions(self):
        return self.get_positions()
    
    @positions.setter
    def positions(self, positions):
        self.set_positions(positions)
    
    def get_positions(self):
        """
        Get global positions.
        """
        from batoms.tools import local2global
        positions = []
        objs = [self.obj, self.obj_o]
        for i in range(2):
            obj = objs[i]
            positions.append(local2global(self.local_positions[i], 
                np.array(obj.matrix_world)))
        return positions
    
    def set_positions(self, vertices, offsets):
        """
        Set global positions to local vertices
        """
        object_mode()
        from batoms.tools import local2global
        n = len(self.obj.data.vertices)
        # scales = np.concatenate((np.ones((n, 1)), 
                        # np.ones((n, 1)),
                        # lengths.reshape(-1, 1)), axis = 1)
        for obj, positions in zip([self.obj, self.obj_o], 
                    [vertices, offsets]):
            if len(positions) != n:
                raise ValueError('positions has wrong shape %s != %s.' %
                                    (len(positions), n))
            positions = local2global(positions, 
                    np.array(obj.matrix_world), reversed = True)
            positions = positions.reshape((n*3, 1))
            obj.data.shape_keys.key_blocks[0].data.foreach_set('co', positions)
            obj.data.update()
        bpy.context.view_layer.objects.active = self.obj
        bpy.ops.object.mode_set(mode = 'EDIT')
        bpy.ops.object.mode_set(mode = 'OBJECT')
    
    @property
    def attributes(self):
        return self.get_attributes()
    
    @attributes.setter
    def attributes(self, attributes):
        self.set_attributes(attributes)

    def get_attributes(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        # attributes
        me = self.obj.data
        nvert = len(me.vertices)
        npoly = len(me.polygons)
        attributes = {}
        for key in me.attributes.keys():
            att = me.attributes.get(key)
            dtype = att.data_type
            domain = att.domain
            n = nvert if domain == 'POINT' else npoly
            if dtype == 'STRING':
                attributes[key] = np.zeros(n, dtype = 'U20')
                for i in range(n):
                    attributes[key][i] = att.data[i].value
            elif dtype == 'INT':
                attributes[key] = np.zeros(n, dtype = int)
                att.data.foreach_get("value", attributes[key])
            elif dtype == 'FLOAT':
                attributes[key] = np.zeros(n, dtype = float)
                att.data.foreach_get("value", attributes[key])
            elif dtype == 'BOOLEAN':
                attributes[key] = np.zeros(n, dtype = bool)
                att.data.foreach_get("value", attributes[key])
            else:
                raise KeyError('%s is not support.'%dtype)
            attributes[key] = np.array(attributes[key])
        return attributes
    
    def set_attributes(self, attributes):
        tstart = time()
        me = self.obj.data
        for key, data in attributes.items():
            # print(key)
            att = me.attributes.get(key)
            if att is None:
                dtype = type(attributes[key][0])
                if np.issubdtype(dtype, int):
                    dtype = 'INT'
                elif np.issubdtype(dtype, float):
                    dtype = 'FLOAT'
                elif np.issubdtype(dtype, str):
                    dtype = 'STRING'
                else:
                    raise KeyError('%s is not supported.'%dtype)
                att = me.attributes.new(name = key, type = dtype, domain = 'POINT')
            if att.data_type == 'STRING':
                nvert = len(me.vertices)
                for i in range(nvert):
                    att.data[i].value = data[i]
            else:
                att.data.foreach_set("value", data)
        me.update()
        # print('set_attributes: %s'%(time() - tstart))

    def set_attribute_with_indices(self, name, indices, data):
        data0 = self.attributes[name]
        data0[indices] = data
        self.set_attributes({name: data0})
    
    
    @property    
    def nframe(self):
        return self.get_nframe()
    
    def get_nframe(self):
        if self.obj.data.shape_keys is None:
            return 0
        nframe = len(self.obj.data.shape_keys.key_blocks)
        return nframe
    
    @property    
    def frames(self):
        return self.get_frames()
    
    @frames.setter    
    def frames(self, frames):
        self.set_frames(frames)
    
    def get_frames(self):
        """
        read shape key
        """
        from batoms.tools import local2global
        obj = self.obj
        n = len(self)
        nframe = self.nframe
        frames = np.empty((nframe, n, 3), dtype=np.float64)
        for i in range(nframe):
            positions = np.empty(n*3, dtype=np.float64)
            sk = obj.data.shape_keys.key_blocks[i]
            sk.data.foreach_get('co', positions)
            local_positions = positions.reshape((n, 3))
            local_positions = local2global(local_positions, 
                            np.array(self.obj.matrix_world))
            frames[i] = local_positions
        return frames
    
    @property    
    def color(self):
        return self.get_color()
    
    @color.setter    
    def color(self, color):
        """
        >>> h.color = [0.8, 0.1, 0.3, 1.0]
        """
        self.set_color(color)
    
    def get_color(self):
        """

        """
        Viewpoint_color = self.material.diffuse_color
        for node in self.material.node_tree.nodes:
            if 'Base Color' in node.inputs:
                node_color = node.inputs['Base Color'].default_value[:]
            if 'Alpha' in node.inputs:
                Alpha = node.inputs['Alpha'].default_value
        color = [node_color[0], node_color[1], node_color[2], Alpha]
        return color
    
    def set_color(self, color):
        if len(color) == 3:
            color = [color[0], color[1], color[2], 1]
        self.material.diffuse_color = color
        for node in self.material.node_tree.nodes:
            if 'Base Color' in node.inputs:
                node.inputs['Base Color'].default_value = color
            if 'Alpha' in node.inputs:
                node.inputs['Alpha'].default_value = color[3]
    
    @property    
    def node(self):
        return self.get_node()
    
    @node.setter    
    def node(self, node):
        self.set_node(node)
    
    def get_node(self):
        return self.material.node_tree.nodes
    
    def set_node(self, node):
        for key, value in node.items():
            self.material.node_tree.nodes['Principled BSDF'].inputs[key].default_value = value
    
    def clean_bpolyhedras_objects(self, obj):
        obj = bpy.data.objects[obj]
        bpy.data.objects.remove(obj, do_unlink = True)
    
    def set_frames(self, frames = None, frame_start = 0, only_basis = False):
        """

        frames: list
            list of positions
        
        >>> from bpolyhedras import Bbond
        >>> import numpy as np
        >>> positions = np.array([[0, 0 ,0], [1.52, 0, 0]])
        >>> h = Bbond('h2o', 'H', positions)
        >>> frames = []
        >>> for i in range(10):
                frames.append(positions + [0, 0, i])
        >>> h.set_frames(frames)
        
        use shape_keys (faster)
        """
        from batoms.butils import add_keyframe_to_shape_key
        from batoms.source_data import bond_source
        from batoms.bdraw import cylinder_mesh_from_vec
        if frames is None:
            frames = self._frames
        vertices, offsets = frames
        nframe = len(vertices)
        if nframe == 0 : return
        nbond = len(vertices[0])
        # scales = np.concatenate((np.ones((nframe, nbond, 1)), 
                        # np.ones((nframe, nbond, 1)),
                        # lengths.reshape(nframe, nbond, 1)), axis = 2)
        for sp, frames in zip(['center', 'offset'],
                                [vertices, offsets]):
            name = '%s_polyhedra_%s'%(self.label, sp)
            obj = bpy.data.objects.get(name)
            base_name = 'Basis_%s_polyhedra_%s'%(self.label, sp)
            if obj.data.shape_keys is None:
                obj.shape_key_add(name = base_name)
            elif base_name not in obj.data.shape_keys.key_blocks:
                obj.shape_key_add(name = base_name)
        if only_basis:
            return
            """
            nvert = len(obj.data.vertices)
            mesh_data = bond_source[self.segments]
            for i in range(1, nframe):
                sk = obj.shape_key_add(name = str(i))
                # Use the local position here
                positions = frames[i][:, 0:3]
                vertices, faces = cylinder_mesh_from_vec(positions[:, 0:3], 
                                    positions[:, 3:6], positions[:, 6:7], 
                                    self.width, mesh_data)
                vertices = vertices.reshape((nvert*3, 1))
                sk.data.foreach_set('co', vertices)
                # Add Keyframes, the last one is different
                if i != nframe - 1:
                    add_keyframe_to_shape_key(sk, 'value', 
                        [0, 1, 0], [frame_start + i - 1, 
                        frame_start + i, frame_start + i + 1])
                else:
                    add_keyframe_to_shape_key(sk, 'value', 
                        [0, 1], [frame_start + i - 1, frame_start + i])

            """
    
    def __len__(self):
        return len(self.obj.data.vertices)
    
    
    def __getitem__(self, index):
        """Return a subset of the Bbond.

        i -- int, describing which atom to return.

        #todo: this is slow for large system
        
        """
        from batoms.bond import Bond
        if isinstance(index, int):
            bond = Bond(self.label, index, polyhedras=self)
            # bpy.ops.object.mode_set(mode=mode)
            return bond
        else:
            return self.positions[index]
        
    
    def __setitem__(self, index, value):
        """Return a subset of the Bbond.

        i -- int, describing which atom to return.

        #todo: this is slow for large system

        """
        positions = self.positions
        positions[index] = value
        self.set_positions(positions)

    
    def repeat(self, m, cell):
        """
        In-place repeat of atoms.

        >>> from batoms.bond import Bbond
        >>> c = Bbond('co', 'C', [[0, 0, 0], [1.2, 0, 0]])
        >>> c.repeat([3, 3, 3], np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]))
        """
        if isinstance(m, int):
            m = (m, m, m)
        for x, vec in zip(m, cell):
            if x != 1 and not vec.any():
                raise ValueError('Cannot repeat along undefined lattice '
                                 'vector')
        M = np.product(m)
        n = len(self)
        positions = np.tile(self.positions, (M,) + (1,) * (len(self.positions.shape) - 1))
        i0 = 0
        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1, m2), cell)
                    i0 = i1
        self.add_vertices(positions[n:])
    
    def copy(self, label, species):
        """
        Return a copy.

        name: str
            The name of the copy.

        For example, copy H species:
        
        >>> h_new = h.copy(label = 'h_new', species = 'H')

        """
        object_mode()
        bbond = Bbond(label, species, self.local_positions, 
                    location = self.obj.location, 
                    scale = self.scale, material=self.material)
        return bbond
    
    def extend(self, other):
        """
        Extend bbond object by appending bbond from *other*.
        
        >>> from batoms.polyhedras import Bbond
        >>> h1 = Bbond('h2o', 'H_1', [[0, 0, 0], [2, 0, 0]])
        >>> h2 = Bbond('h2o', 'H_2', [[0, 0, 2], [2, 0, 2]])
        >>> h = h1 + h2
        """
        # could also use self.add_vertices(other.positions)
        object_mode()
        bpy.ops.object.select_all(action='DESELECT')
        self.obj.select_set(True)
        other.obj.select_set(True)
        bpy.context.view_layer.objects.active = self.obj
        bpy.ops.object.join()
    
    def __iadd__(self, other):
        """
        >>> h1 += h2
        """
        self.extend(other)
        return self
    
    def __add__(self, other):
        """
        >>> h1 = h1 + h2
        """
        self += other
        return self
    
    def __iter__(self):
        bbond = self.obj
        for i in range(len(self)):
            yield bbond.matrix_world @ bbond.data.vertices[i].co
    
    def __repr__(self):
        s = "polyhedras(Total: {:6d}, {}" .format(len(self), self.arrays)
        return s
    
    def add_vertices(self, positions):
        """
        Todo: find a fast way.
        """
        object_mode()
        positions = positions - self.location
        bm = bmesh.new()
        bm.from_mesh(self.obj.data)
        bm.verts.ensure_lookup_table()
        verts = []
        for pos in positions:
            bm.verts.new(pos)
        bm.to_mesh(self.obj.data)

    def calc_polyhedra_data(self, bondlists, positions, cell):
        """
        """
        from ase.data import chemical_symbols
        from batoms.tools import get_polyhedra_kind
        from scipy.spatial import ConvexHull, Delaunay

        tstart = time()
        chemical_symbols = np.array(chemical_symbols)
        # total 
        npoly = 0
        for b in self.setting:
            sp = string2Number(b.species)
            ind = np.where(bondlists['species_index1'] == sp)[0]
            npoly += len(ind)
        if npoly == 0:
            return default_polyhedra_datas
        # properties
        npoly = npoly*12
        atoms_index1 = np.zeros(npoly, dtype = int) # 0, 1, 2
        atoms_index2 = np.zeros(npoly, dtype = int) # 0, 1, 2
        face_species_index = np.zeros(npoly, dtype = int) # 0, 1, 2
        styles = np.zeros(npoly, dtype = int) # 0, 1, 2
        widths = np.ones(npoly, dtype = float) 
        species_index = np.zeros(npoly, dtype = int)
        #------------------------------------
        offsets0 = bondlists['offsets']
        positions0 = positions[bondlists['atoms_index2']] + np.dot(offsets0, cell)
        #---------------------------------------------
        vertices = []
        offsets = []
        edges = []
        faces = []
        nv = 0
        nf = 0
        for poly in self.setting:
            sp = string2Number(poly.species)
            # loop center atoms
            spis = np.where(bondlists['species_index1'] == sp)[0]
            indices1 = bondlists['atoms_index1'][spis]
            argsort = indices1.argsort()
            indices1 = indices1[argsort]
            indices2 = bondlists['atoms_index2'][spis][argsort]
            positions1 = positions0[spis][argsort]
            offsets1 = offsets0[spis][argsort]
            u, indices = np.unique(indices1, return_index=True)
            indices = np.append(indices, len(indices1))
            n = len(u)
            # print(n)
            for i in range(n):
                vertices1 = positions1[indices[i]: indices[i + 1]]
                dnv = len(vertices1)
                if dnv >= 4:
                    # search face for polyhedra
                    hull = ConvexHull(vertices1)
                    face = hull.simplices
                    # tri = Delaunay(vertices1)
                    # face = tri.simplices
                    dnf = len(face)
                    face = face + nv
                    edge = []
                    for f in face:
                        edge.append([f[0], f[1]])
                        edge.append([f[0], f[2]])
                        edge.append([f[1], f[2]])
                    edges = edges + list(edge)
                    faces = faces + list(face)
                    vertices.extend(vertices1)
                    offsets.extend(offsets1[indices[i]: indices[i + 1]])
                    nv1 = nv + dnv
                    nf1 = nf + dnf
                    atoms_index1[nv:nv1] = indices1[indices[i]: indices[i + 1]]
                    atoms_index2[nv:nv1] = indices2[indices[i]: indices[i + 1]]
                    styles[nv:nv1] = int(poly.style)
                    widths[nv:nv1] = poly.width
                    species_index[nv:nv1] = string2Number(poly.species)
                    face_species_index[nf:nf1] = string2Number(poly.species)
                    nv = nv1
                    nf = nf1
        n = len(vertices)
        nf = len(faces)
        datas = {
            'atoms_index1': atoms_index1[0:n],
            'atoms_index2': atoms_index2[0:n],
            'species_index': species_index[0:n],
            'face_species_index': face_species_index[0:nf],
            'vertices':np.array(vertices),
            'offsets':np.array(offsets),
            'widths':widths[0:n],
            'edges':edges,
            'faces':faces,
            'styles':styles[0:n],
        }
        print('datas: ', datas)
        print('calc_polyhedra_data: {0:10.2f} s'.format(time() - tstart))
        return datas
