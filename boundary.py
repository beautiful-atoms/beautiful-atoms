from pandas import array
import bpy
from ase import Atoms
import numpy as np
from time import time
from batoms.base import ObjectGN
# from batoms.tools import build_boundary
from batoms.butils import object_mode
from batoms.tools import number2String, string2Number
from ase.geometry import wrap_positions, complete_cell


shapes = ["UV_SPHERE", "ICO_SPHERE", "CUBE", "METABALL"]

default_attributes = [
        ['atoms_index', 'INT'], 
        ['species_index', 'INT'], 
        ['show', 'BOOLEAN'], 
        ['select', 'INT'],
        ['model_style', 'INT'],
        ['scale', 'FLOAT'], 
        ['radius_style', 'INT'],
        ]

default_boundary_datas = {
        'atoms_index': np.ones(0, dtype = int),
        'species_index': np.ones(0, dtype = int),
        'species': np.ones(0, dtype = 'U4'),
        'positions':np.zeros((0, 3)),
        'scales':np.zeros(0),
        'offsets':np.zeros((0, 3)),
        'model_styles':np.ones(0, dtype = int),
        'radius_styles':np.ones(0, dtype = int),
        'shows':np.ones(0, dtype = int),
        'selects':np.ones(0, dtype = int),
        }

boundary0 = np.array([[-0.01, 1.01], [-0.01, 1.01], [-0.01, 1.01]])


def build_boundary(species, positions, cell, pbc, boundary, include_self = False):
    """
    find region outside cell for periodic caculation, supercell 
    
    todo: boundary > 1
    """
    from functools import reduce
    
    tstart = time()
    if isinstance(boundary, float):
        boundary = [[-boundary, 1 + boundary], [-boundary, 1+boundary], [-boundary, 1+boundary]]
    boundary = np.array(boundary)
    ib = np.concatenate((np.ceil(boundary[:, 0]).reshape(-1, 1), 
                np.floor(boundary[:, 1]).reshape(-1, 1)), axis = 1)
    delta = boundary - ib
    wraped_positions = wrap_positions(positions, cell, pbc=pbc)
    scaled_positions = np.linalg.solve(complete_cell(cell).T,
                                              wraped_positions.T).T
    # boundary condition
    natom = len(positions)
    ind = np.arange(natom)
    # find atoms close to cell face with distance of r
    indices = [[], [], []]
    nb = 0
    for i in range(3):
        if not pbc[i]: continue
        ind0 =  np.where((scaled_positions[:, i] < delta[i, 1]))[0]
        ind1 =  np.where((scaled_positions[:, i] > 1- delta[i, 0]))[0]
        indices[i] = [ind0, ind1]
        nb += len(ind0) + len(ind1)
    print(indices)
    
def search_boundary(species, positions, cell, pbc, boundary = [[0, 1], [0, 1], [0, 1]], skin = 3):
    """
    search atoms in the boundary

    Parameters:

    atoms: ASE Atoms object

    boundary: list

    skin: float
        Could be the maximum cutoff.

    Return:

    atoms_boundary: ASE Atoms object
        atoms inside the boudary, not include the core
    offsets_skin: numpy array
        atoms close to boundary with a skin distance. 
        Used to search bond outside the boundary.
        index and offset

    """
    tstart = time()
    if isinstance(boundary, float):
        boundary = [[-boundary, 1 + boundary], [-boundary, 1+boundary], [-boundary, 1+boundary]]
    boundary = np.array(boundary)
    boundary_skin = boundary.copy()
    # skin to scaled distance in cell
    par = np.linalg.norm(cell, axis = 0)
    skin = np.array([skin/par[0], skin/par[1], skin/par[2]])
    # skin region
    boundary_skin[:, 0] += skin
    boundary_skin[:, 1] -= skin
    # find supercell
    f = np.floor(boundary)
    c = np.ceil(boundary)
    ib = np.array([f[:, 0], c[:, 1]]).astype(int)
    M = np.product(ib[1] - ib[0] + 1)
    positions = np.linalg.solve(complete_cell(cell).T,
                                              positions.T).T
    n = len(positions)
    npositions = np.tile(positions, (M - 1,) + (1,) * (len(positions.shape) - 1))
    i0 = 0
    # index
    offsets = np.zeros((M*n, 4), dtype=int)
    ind0 = np.arange(n).reshape(-1, 1)
    # indices = np.tile(ind0, (M - 1,) + (1,) * (len(ind0.shape) - 1))
    # symbols0 = elements
    species0 = species
    symbols = []
    species = []
    # repeat the positions so that
    # it completely covers the boundary
    for m0 in range(ib[0, 0], ib[1, 0] + 1):
        for m1 in range(ib[0, 1], ib[1, 1] + 1):
            for m2 in range(ib[0, 2], ib[1, 2] + 1):
                if m0 == 0 and m1 == 0 and m2 == 0: continue
                i1 = i0 + n
                npositions[i0:i1] += (m0, m1, m2)
                offsets[i0:i1] = np.append(ind0, np.array([[m0, m1, m2]]*n), axis = 1)
                # symbols.extend(symbols0)
                species.extend(species0)
                i0 = i1
    # boundary condition
    ind1 =  np.where((npositions[:, 0] > boundary[0][0]) & (npositions[:, 0] < boundary[0][1]) \
                & (npositions[:, 1] > boundary[1][0]) & (npositions[:, 1] < boundary[1][1]) \
                & (npositions[:, 2] > boundary[2][0]) & (npositions[:, 2] < boundary[2][1]))[0]
    # build atoms inside the boundary
    # indices_b = indices[ind1]
    # npositions_b = npositions[ind1]
    # npositions_b = np.dot(npositions_b, cell)
    # symbols_b = np.array(symbols)[ind1]
    # species_b = np.array(species)[ind1]
    offsets_b = offsets[ind1]
    # print('search boundary: {0:10.2f} s'.format(time() - tstart))
    return offsets_b

def search_bond(positions0, offsets_skin, bondlists, boundary, recursive = False, previous = None):
    """
    Search bonded atoms of sp1 or sp2 recursively.

    Parameters:

    positions0: array
        positions of original atoms, the core atoms
    offsets_skin: array
        index and offset of atoms of sp1 or sp2
    bondlists: array
        bondlist of the core atoms
    boundary: array
    recursive: bool
        recursive or not, for search mode 1 or 2
    previous: array
        offsets_skin of previous search. 
        The index should be remove from the result of this search.
    
    Return:

    offset_new: array
        index and offset of atoms which bond to sp1 or sp2
    
    """
    neighbors = np.array([]).reshape(-1, 4)
    if len(bondlists) ==0 or len(offsets_skin) == 0: return neighbors
    # build bonded atoms using offset_skin and bondlists
    sites = set(offsets_skin[:, 0])
    for i in sites:
        ind = np.where(bondlists[:, 0] == i)[0]
        if len(ind) == 0: continue
        neighbor = bondlists[ind][:, 1:]
        sitei = np.where(offsets_skin[:, 0] == i)[0]
        for j in sitei:
            temp = neighbor.copy()
            temp[:, 1:4] = temp[:, 1:4] + offsets_skin[j][1:4]
            neighbors = np.append(neighbors, temp, axis = 0)
    # choose index by two pricinple
    # 1, outside the boundary
    # 2, not from preious index
    # rebuild the positions, and check boundary
    neighbors = neighbors.astype(int)
    npositions = positions0[neighbors[:, 0]] + neighbors[:, 1:4]
    index =  np.where((npositions[:, 0] > boundary[0][0]) & (npositions[:, 0] < boundary[0][1]) \
                & (npositions[:, 1] > boundary[1][0]) & (npositions[:, 1] < boundary[1][1]) \
                & (npositions[:, 2] > boundary[2][0]) & (npositions[:, 2] < boundary[2][1]))
    #
    mask = np.ones(npositions.shape[0], dtype=bool)
    mask[index] = False
    offset_new = neighbors[mask]
    # offset_new = np.unique(offset_new, axis=0)
    # remove previous index, and remove duplicate
    if previous is not None:
        __, indices = np.unique(np.concatenate([previous, offset_new]), return_index=True, axis=0)
        indices = indices[indices >= len(previous)] - len(previous)
        offset_new = offset_new[indices]
    if len(offset_new) == 0: return np.array([]).reshape(-1, 4)
    # recursive search
    if recursive:
        offset_new1 = search_bond(positions0, offset_new, bondlists, boundary, recursive = True, previous = offsets_skin)
        # save
        if offset_new1 is not None:
            offset_new = np.append(offset_new, offset_new1, axis = 0)
    else:
        return offset_new
    return offset_new


class Boundary(ObjectGN):
    """Boundary Class

    """
    def __init__(self, 
                label = None,
                boundary = np.array([[-0.0, 1.0], [-0.0, 1.0], [-0.0, 1.0]]),
                boundary_datas = None,
                batoms = None,
                 ):
        #
        self.batoms = batoms
        self.label = label
        name = 'boundary'
        ObjectGN.__init__(self, label, name)
        if boundary_datas is not None:
            self.build_object(boundary_datas)
        else:
            self.load(label)
        self.boundary = boundary

    def build_object(self, boundary_datas, location = (0, 0, 0), attributes = {}):
        """
        build child object and add it to main objects.
        """
        tstart = time()
        if len(boundary_datas['positions'].shape) == 2:
            self._frames = {'positions': np.array([boundary_datas['positions']]),
                            'offsets': np.array([boundary_datas['offsets']]),
                        }
            positions = boundary_datas['positions']
            offsets = boundary_datas['offsets']
        elif len(boundary_datas['positions'].shape) == 3:
            self._frames = {'positions': boundary_datas['positions'],
                            'offsets': boundary_datas['offsets'],
                            }
            positions = boundary_datas['positions'][0]
            offsets = boundary_datas['offsets'][0]
        else:
            raise Exception('Shape of positions is wrong!')
        #
        attributes.update({
                            'atoms_index': boundary_datas['atoms_index'], 
                            'species_index': boundary_datas['species_index'], 
                            'show': boundary_datas['shows'],
                            'model_style': boundary_datas['model_styles'],
                            'select': boundary_datas['selects'],
                            'scale': boundary_datas['scales'],
                            'radius_style': boundary_datas['radius_styles'],
                            })
        name = '%s_boundary'%self.label
        self.delete_obj(name)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(positions, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        for attribute in default_attributes:
            mesh.attributes.new(name = attribute[0], type = attribute[1], domain = 'POINT')
        self.batoms.coll.objects.link(obj)
        obj.parent = self.batoms.obj
        #
        name = '%s_boundary_offset'%self.label
        self.delete_obj(name)
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(offsets, [], [])
        mesh.update()
        obj = bpy.data.objects.new(name, mesh)
        self.batoms.coll.objects.link(obj)
        obj.hide_set(True)
        obj.parent = self.batoms.obj
        bpy.context.view_layer.update()
        self.set_attributes(attributes)
        self.build_geometry_node()
        self.set_frames(self._frames, only_basis = True)
        # print('boundary: build_object: {0:10.2f} s'.format(time() - tstart))
    
    def build_geometry_node(self):
        """
        """
        from batoms.butils import get_nodes_by_name
        name = 'GeometryNodes_%s_boundary'%self.label
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
        for i in range(1, 7):
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
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                        '%s_JoinGeometry'%self.label, 
                        'GeometryNodeJoinGeometry')
        gn.node_group.links.new(GroupInput.outputs['Geometry'], JoinGeometry.inputs['Geometry'])
        gn.node_group.links.new(JoinGeometry.outputs['Geometry'], GroupOutput.inputs['Geometry'])
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
        gn.node_group.links.new(GroupInput.outputs[1], TransferBatoms.inputs['Index'])
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
        OffsetNode = self.vectorDotMatrix(gn, TransferOffsets, self.batoms.cell, '')
        # we need one add operation to get the positions with offset
        VectorAdd = get_nodes_by_name(gn.node_group.nodes, 
                    '%s_VectorAdd'%(self.label),
                    'ShaderNodeVectorMath')
        VectorAdd.operation = 'ADD'
        gn.node_group.links.new(TransferBatoms.outputs[0], VectorAdd.inputs[0])
        gn.node_group.links.new(OffsetNode.outputs[0], VectorAdd.inputs[1])
        # set positions
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                        '%s_SetPosition'%self.label, 
                        'GeometryNodeSetPosition')
        gn.node_group.links.new(GroupInput.outputs['Geometry'], SetPosition.inputs['Geometry'])
        gn.node_group.links.new(VectorAdd.outputs[0], SetPosition.inputs['Position'])
        
    def add_geometry_node(self, spname, selname):
        """
        """
        from batoms.butils import get_nodes_by_name
        gn = self.gnodes
        GroupInput = gn.node_group.nodes.get('Group Input')
        SetPosition = get_nodes_by_name(gn.node_group.nodes,
                        '%s_SetPosition'%self.label)
        JoinGeometry = get_nodes_by_name(gn.node_group.nodes,
                        '%s_JoinGeometry'%self.label, 
                        'GeometryNodeJoinGeometry')
        CompareSelect = get_nodes_by_name(gn.node_group.nodes, 
                    'select_%s_%s'%(self.label, selname),
                    'FunctionNodeCompareFloats')
        CompareSelect.operation = 'EQUAL'
        # CompareSelect.data_type = 'INT'
        CompareSelect.inputs[1].default_value = string2Number(selname)
        gn.node_group.links.new(GroupInput.outputs[4], CompareSelect.inputs[0])
        CompareSpecies = get_nodes_by_name(gn.node_group.nodes, 
                    'CompareFloats_%s_%s'%(self.label, spname),
                    'FunctionNodeCompareFloats')
        CompareSpecies.operation = 'EQUAL'
        # CompareSpecies.data_type = 'INT'
        CompareSpecies.inputs[1].default_value = string2Number(spname)
        InstanceOnPoint = get_nodes_by_name(gn.node_group.nodes,
                    'InstanceOnPoint_%s_%s_%s'%(self.label, selname, spname), 
                    'GeometryNodeInstanceOnPoints')
        ObjectInfo = get_nodes_by_name(gn.node_group.nodes, 
                    'ObjectInfo_%s_%s_%s'%(self.label, selname, spname),
                    'GeometryNodeObjectInfo')
        ObjectInfo.inputs['Object'].default_value = self.batoms.species.instancers[selname][spname]
        #
        BoolSelectSpecies = get_nodes_by_name(gn.node_group.nodes, 
                        'BooleanMath_%s_%s_%s_0'%(self.label, selname, spname),
                        'FunctionNodeBooleanMath')
        BoolShow = get_nodes_by_name(gn.node_group.nodes, 
                    'BooleanMath_%s_%s_%s_1'%(self.label, selname, spname),
                    'FunctionNodeBooleanMath')
        #
        gn.node_group.links.new(SetPosition.outputs['Geometry'], InstanceOnPoint.inputs['Points'])
        gn.node_group.links.new(GroupInput.outputs[2], CompareSpecies.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[3], BoolShow.inputs[0])
        gn.node_group.links.new(GroupInput.outputs[6], InstanceOnPoint.inputs['Scale'])
        gn.node_group.links.new(CompareSelect.outputs[0], BoolSelectSpecies.inputs[0])
        gn.node_group.links.new(CompareSpecies.outputs[0], BoolSelectSpecies.inputs[1])
        gn.node_group.links.new(BoolSelectSpecies.outputs[0], BoolShow.inputs[1])
        gn.node_group.links.new(BoolShow.outputs['Boolean'], InstanceOnPoint.inputs['Selection'])
        gn.node_group.links.new(ObjectInfo.outputs['Geometry'], InstanceOnPoint.inputs['Instance'])
        gn.node_group.links.new(InstanceOnPoint.outputs['Instances'], JoinGeometry.inputs['Geometry'])        
    
    def update(self):
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
        tstart = time()
        for f in range(nframe):
            # print('update boundary: ', f)
            positions = frames[f]
            boundary_list = search_boundary(species, positions, self.batoms.cell,
            self.batoms.pbc, self.boundary)
            # search bond
            bondlist = self.bondlists
            """
            if len(bondlist['positions']) == 0:
                self.batoms.bonds.update()
                bondlist = self.batoms.bonds.arrays
            #
            offsets_skin1, bondlist1, offsets_skin2, bondlist2 = \
                    self.bondsetting.search_bond_list(atoms0, bondlists, offsets_skin)
            # search type 1
            offset_skin_1 = search_bond(atoms0.get_scaled_positions(), offsets_skin1, bondlist1, boundary)
            # search type 2
            offset_skin_2 = search_bond(atoms0.get_scaled_positions(), offsets_skin2, bondlist2, boundary, recursive=True)
            # search type 1 final
            offset_skin_3 = search_bond(atoms0.get_scaled_positions(), offset_skin_2, bondlist1, boundary)
            """
            # print('calc bond: {0:10.2f} s'.format(time() - tstart))
            
            if f == 0:
                boundary_lists = boundary_list
            else:
                boundary_lists = np.append(boundary_lists, boundary_list, axis = 0)
        boundary_datas = self.calc_boundary_data(boundary_lists, arrays, self.batoms.cell)
        # if len(boundary_datas) == 0:
            # return
        self.set_arrays(boundary_datas)
        self.batoms.draw()
        # self.coll.objects.link(bb.obj)
        # bpy.data.collections['Collection'].objects.unlink(bb.obj)
        # bb.set_frames()
        # bpy.context.scene.frame_set(self.batoms.nframe)
        print('draw bond: {0:10.2f} s'.format(time() - tstart))
    
    @property
    def obj_o(self):
        return self.get_obj_o()
    
    def get_obj_o(self):
        name = '%s_boundary_offset'%self.label
        obj_o = bpy.data.objects.get(name)
        if obj_o is None:
            raise KeyError('%s object is not exist.'%name)
        return obj_o
    
    @property    
    def boundary(self):
        return self.get_boundary()
    
    @boundary.setter    
    def boundary(self, boundary):
        if boundary is not None:
            if isinstance(boundary, (int, float)):
                boundary = np.array([[-boundary, 1 + boundary]]*3)
            elif len(boundary) == 3:
                if isinstance(boundary[0], (int, float)):
                    boundary = np.array([[-boundary[0], 1 + boundary[0]],
                                      [-boundary[1], 1 + boundary[1]],
                                      [-boundary[2], 1 + boundary[2]]])
                elif len(boundary[0]) == 2:
                    boundary = np.array(boundary)
            else:
                raise Exception('Wrong boundary setting!')
            self.batoms.coll.batoms.boundary = boundary[:].flatten()
        self.update()
        # if self.model_type == 1:
        #     self.draw_bonds()
        # if self.model_type == 2:
        #     self.draw_bonds()
        #     self.draw_polyhedras()
    
    def get_boundary(self):
        boundary = np.array(self.batoms.coll.batoms.boundary)
        return boundary.reshape(3, -1)

    def set_arrays(self, arrays):
        """
        """
        # if len(arrays['positions']) == 0:
            # return
        attributes = self.attributes
        # same length
        if len(arrays['positions']) == len(attributes['show']):
            self.positions = arrays['positions']
            self.offsets = arrays['offsets']
            species_index = [string2Number(sp) for sp in arrays['species']]
            self.set_attributes({'species_index': species_index,
                                'scale': arrays['scales'],
                                'show': arrays['shows'],
                                })
        else:
            # add or remove vertices
            self.build_object(arrays)
            species = np.unique(arrays['species'])
            for sp in species:
                self.add_geometry_node(sp, 'sel0')

    def get_arrays(self):
        """
        """
        object_mode()
        tstart = time()
        arrays = self.attributes
        arrays.update({'positions': self.positions})
        arrays.update({'offsets': self.offsets})
        # radius
        radius = self.batoms.radius
        arrays.update({'radius': np.zeros(len(self))})
        species = np.array([number2String(i) for i in arrays['species_index']], dtype = 'U20')
        arrays['species'] = species
        for sel, data in radius.items():
            for sp, value in data.items():
                mask = np.where((arrays['species'] == sp) & (arrays['select'] == string2Number(sel)))
                arrays['radius'][mask] = value
        # size
        arrays['size'] = arrays['radius']*arrays['scale']
        # main elements
        main_elements = self.batoms.species.main_elements
        elements = [main_elements[sp] for sp in arrays['species']]
        arrays.update({'elements': np.array(elements, dtype='U20')})
        # print('get_arrays: %s'%(time() - tstart))
        return arrays

    @property
    def boundary_data(self):
        return self.get_boundary_data()

    def get_boundary_data(self, include_batoms = False):
        """
        using foreach_get and foreach_set to improve performance.
        """
        arrays = self.arrays
        boundary_data = {'positions': arrays['positions'],
                        'species': arrays['species'],
                        'indices': arrays['atoms_index'],
                        'offsets': arrays['offsets']}
        return boundary_data

    @property
    def offsets(self):
        return self.get_offsets()
    
    def get_offsets(self):
        """
        using foreach_get and foreach_set to improve performance.
        """
        n = len(self)
        offsets = np.empty(n*3, dtype = int)
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
        bpy.context.view_layer.objects.active = self.obj_o
        bpy.ops.object.mode_set(mode = 'EDIT')
        bpy.ops.object.mode_set(mode = 'OBJECT')
    
    @property
    def bondlists(self):
        return self.get_bondlists()
    
    def get_bondlists(self):
        try:
            bondlists = self.batoms.bonds.arrays
        except:
            bondlists = None
        return bondlists

    def get_frames(self):
        """
        """
        frames = {}
        frames['positions'] = self.get_obj_frames(self.obj)
        frames['offsets'] = self.get_obj_frames(self.obj_o)
        return frames

    def set_frames(self, frames = None, frame_start = 0, only_basis = False):
        if frames is None:
            frames = self._frames
        nframe = len(frames)
        if nframe == 0 : return
        name = '%s_boundary'%(self.label)
        obj = self.obj
        self.set_obj_frames(name, obj, frames['positions'])
        #
        name = '%s_boundary_offset'%(self.label)
        obj = self.obj_o
        self.set_obj_frames(name, obj, frames['offsets'])

    def calc_boundary_data(self, boundary_lists, arrays, cell):
        """
        """

        tstart = time()
        # properties
        nb =len(boundary_lists)
        model_styles = arrays['model_style'][boundary_lists[:, 0]]
        shows = arrays['show'][boundary_lists[:, 0]]
        radius_styles = arrays['radius_style'][boundary_lists[:, 0]]
        selects = arrays['select'][boundary_lists[:, 0]]
        scales = arrays['scale'][boundary_lists[:, 0]]
        species_indexs = arrays['species_index'][boundary_lists[:, 0]]
        species = arrays['species'][boundary_lists[:, 0]]
        #------------------------------------
        offset_vectors = boundary_lists[:, 1:4]
        offsets = np.dot(offset_vectors, cell)
        positions = arrays['positions'][boundary_lists[:, 0]] + offsets
        datas = {
            'atoms_index': np.array(boundary_lists[:, 0]),
            'species_index': species_indexs,
            'species': species,
            'positions':positions,
            # 'offsets':offsets,
            'offsets':offset_vectors,
            'model_styles':model_styles,
            'shows':shows,
            'selects':selects,
            'scales':scales,
            'radius_styles':radius_styles,
        }
        # print('datas: ', datas)
        # print('calc_boundary_data: {0:10.2f} s'.format(time() - tstart))
        return datas

    def __getitem__(self, index):
        return self.boundary[index]

    def __setitem__(self, index, value):
        """Set boundary vectors.
        """
        boundary = self.boundary
        if isinstance(value, (int, float)):
            boundary = np.array([-boundary, 1 + boundary])
        boundary[index] = value
        self.boundary = boundary
    def __repr__(self) -> str:
        s = self.boundary.__repr__()
        return s