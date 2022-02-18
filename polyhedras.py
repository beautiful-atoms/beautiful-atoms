"""Definition of the polyhedras class.

This module defines the polyhedras object in the Batoms package.

"""

import bpy
import bmesh
from time import time
from batoms.butils import object_mode, compareNodeType
from batoms.tools import string2Number
import numpy as np
from batoms.base import ObjectGN
from batoms.polyhedrasetting import PolyhedraSettings

default_attributes = [
            ['atoms_index1', 'INT', 'POINT'],
            ['atoms_index2', 'INT', 'POINT'],
            ['species_index', 'INT', 'POINT'],
            ['face_species_index', 'INT', 'FACE'],
            ['show', 'BOOLEAN', 'POINT'],
            ['style', 'INT', 'POINT'],
        ]

default_GroupInput = [
            ['atoms_index1', 'NodeSocketInt', 'POINT'],
            ['atoms_index2', 'NodeSocketInt', 'POINT'],
            ['species_index', 'NodeSocketInt', 'POINT'],
            ['face_species_index', 'NodeSocketInt', 'FACE'],
            ['show', 'NodeSocketBool', 'POINT'],
            ['style', 'NodeSocketInt', 'POINT'],
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

class Polyhedras(ObjectGN):
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
        name = 'polyhedra'
        ObjectGN.__init__(self, label, name)
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
            self._frames = {'vertices': np.array([polyhedra_datas['vertices']]),
                            'offsets': np.array([polyhedra_datas['offsets']]),
                        }
            vertices = polyhedra_datas['vertices']
            offsets = polyhedra_datas['offsets']
        elif len(polyhedra_datas['vertices'].shape) == 3:
            self._frames = {'vertices': polyhedra_datas['vertices'],
                            'offsets': polyhedra_datas['offsets'],
                            }
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
        name = self.obj_name
        self.delete_obj(name)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(vertices, polyhedra_datas['edges'], polyhedra_datas['faces'])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        for attribute in default_attributes:
            mesh.attributes.new(name = attribute[0], type = attribute[1], domain = attribute[2])
        self.setting.coll.objects.link(obj)
        obj.parent = self.batoms.obj
        #
        name = '%s_polyhedra_offset'%self.label
        self.delete_obj(name)
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
        obj.parent = self.batoms.obj
        print('polyhedras: build_object: {0:10.2f} s'.format(time() - tstart))

    def assign_materials(self):
        # sort element by occu
        me = self.obj.data
        me.materials.clear()
        for sp in self.setting:
            me.materials.append(self.setting.materials[sp.species][0])
    
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
        inputs = modifier.node_group.inputs
        GroupInput = modifier.node_group.nodes.get('Group Input')
        GroupOutput = modifier.node_group.nodes.get('Group Output')
        # add new output sockets
        for att in default_GroupInput:
            GroupInput.outputs.new(type = att[1], name = att[0])
            inputs.new(att[1], att[0])
            id = inputs[att[0]].identifier
            modifier['%s_use_attribute'%id] = True
            modifier['%s_attribute_name'%id] = att[0]
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
        gn.node_group.links.new(ObjectBatoms.outputs['Geometry'], TransferBatoms.inputs[0])
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
        gn.node_group.links.new(ObjectOffsets.outputs['Geometry'], TransferOffsets.inputs[0])
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
                    compareNodeType)
        CompareSpecies.operation = 'EQUAL'
        CompareSpecies.inputs[1].default_value = string2Number(sp["species"])
        gn.node_group.links.new(GroupInput.outputs[3], CompareSpecies.inputs[0])
        # face
        CompareSpeciesFace = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_CompareSpecies_%s_face'%(self.label, sp["species"]),
                    compareNodeType)
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

    def update(self, ):
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
        show = arrays['show']
        species = arrays['species']
        # frames_boundary = self.batoms.get_frames(self.batoms.batoms_boundary)
        # frames_search = self.batoms.get_frames(self.batoms.batoms_search)
        nframe = len(frames)
        polyhedra_datas = {}
        tstart = time()
        for f in range(nframe):
            print('update polyhedra: ', f)
            positions = frames[f, show, :]
            # if len(frames_boundary) > 0:
            #     positions_boundary = frames_boundary[f]
            #     positions = positions + positions_boundary
            # if len(frames_search) > 0:
            #     positions_search = frames_search[f]
            #     positions = positions + positions_search
            bondlists = self.bondlists
            if len(bondlists) == 0:
                self.batoms.bonds.update()
                bondlists = self.batoms.bonds.bondlists
            polyhedra_kinds = self.calc_polyhedra_data(bondlists, species, 
                        positions, arrays['model_style'][show])
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
        self.set_arrays(polyhedra_datas)
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
            bondlists = self.batoms.bonds.bondlists
        except:
            bondlists = None
        return bondlists
        
    def set_arrays(self, arrays):
        """
        """
        attributes = self.attributes
        # same length
        if len(arrays['vertices']) == len(attributes['show']):
            self.positions = arrays['vertices']
            self.offsets = arrays['offsets']
            self.set_attributes({'atoms_index1': arrays['atoms_index1']})
            self.set_attributes({'atoms_index2': arrays['atoms_index2']})
            self.set_attributes({'species_index': arrays['species_index']})
            self.set_attributes({'face_species_index': arrays['face_species_index']})
            self.set_attributes({'style': arrays['styles']})
        else:
            # add or remove vertices, rebuild the object
            self.build_object(arrays)

    @property
    def obj(self):
        return self.get_obj()
    
    def get_obj(self):
        obj = bpy.data.objects.get(self.obj_name)
        if obj is None:
            self.build_object(default_polyhedra_datas)
        return obj

    @property
    def offsets(self):
        return self.get_offsets()
    
    def get_offsets(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        n = len(self)
        offsets = np.empty(n*3, dtype=np.float64)
        self.obj_o.data.vertices.foreach_get('co', offsets)  
        return offsets.reshape((n, 3))
    
    @offsets.setter
    def offsets(self, offsets):
        self.set_offsets(offsets)
    
    def set_offsets(self, offsets):
        """
        Set global offsets to local vertices
        """
        object_mode()
        from batoms.tools import local2global
        n = len(self.obj_o.data.vertices)
        if len(offsets) != n:
            raise ValueError('offsets has wrong shape %s != %s.' %
                                (len(offsets), n))
        offsets = offsets.reshape((n*3, 1))
        self.obj_o.data.shape_keys.key_blocks[0].data.foreach_set('co', offsets)
        self.obj_o.data.update()
        
    
    def get_frames(self):
        """
        """
        frames = {}
        frames['vertices'] = self.get_obj_frames(self.obj)
        frames['offsets'] = self.get_obj_frames(self.obj_o)
        return frames
        
    def set_frames(self, frames = None, frame_start = 0, only_basis = False):
        if frames is None:
            frames = self._frames
        nframe = len(frames)
        if nframe == 0 : return
        name = '%s_polyhedra'%(self.label)
        obj = self.obj
        self.set_obj_frames(name, obj, frames['vertices'])
        #
        name = '%s_polyhedra_offset'%(self.label)
        obj = self.obj_o
        self.set_obj_frames(name, obj, frames['offsets'])
        
    
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

    def calc_polyhedra_data(self, bondlists, species, positions,
                        model_styles):
        """
        """
        from ase.data import chemical_symbols
        from batoms.tools import get_polyhedra_kind
        from scipy.spatial import ConvexHull, Delaunay

        tstart = time()
        chemical_symbols = np.array(chemical_symbols)
        # find bonds contribute to polyhedra
        indices = bondlists[:, 8].astype(bool)
        bondlists = bondlists[indices]
        # maxinum number of poly 
        npoly = len(bondlists)
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
        offsets10 = bondlists[:, 2:5]
        offsets20 = bondlists[:, 5:8]
        indices10 = bondlists[:, 0].astype(int)
        indices20 = bondlists[:, 1].astype(int)
        positions10 = positions[indices10] + offsets10
        positions20 = positions[indices20] + offsets20
        #---------------------------------------------
        vertices = []
        offsets = []
        edges = []
        faces = []
        nv = 0
        nf = 0
        speciesarray = species[indices10]
        for poly in self.setting:
            # find center atoms == species
            spis = np.where(speciesarray == poly.species)[0]
            indices1 = indices10[spis]
            bondlist1 = bondlists[spis]
            # center atoms is define by i and the offset
            atom_centers = bondlist1[:, [0, 2, 3, 4]]
            positions21 = positions20[spis]
            # find indices of center atoms
            u, indices = np.unique(atom_centers, axis = 0, return_inverse = True)
            indices = np.append(indices, len(indices1))
            n = len(u)
            # print(n)
            for i in range(n):
                if model_styles[int(u[i, 0])] != 2: continue
                mask = np.where(indices == i)[0]
                vertices1 = positions21[mask]
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
                    offsets.extend(bondlist1[mask, 5:8])
                    nv1 = nv + dnv
                    nf1 = nf + dnf
                    atoms_index1[nv:nv1] = bondlist1[mask, 0]
                    atoms_index2[nv:nv1] = bondlist1[mask, 1]
                    styles[nv:nv1] = int(poly.style)
                    widths[nv:nv1] = poly.width
                    species_index[nv:nv1] = string2Number(poly.species)
                    face_species_index[nf:nf1] = string2Number(poly.species)
                    nv = nv1
                    nf = nf1
        n = len(vertices)
        if n == 0:
            datas = default_polyhedra_datas
        else:
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
                'model_styles':model_styles,
            }
        # print('datas: ', datas)
        print('calc_polyhedra_data: {0:10.2f} s'.format(time() - tstart))
        return datas
