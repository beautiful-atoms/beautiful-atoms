"""
"""
import bpy
from batoms.butils import object_mode
import numpy as np
from time import time
from batoms.base import Setting, tuple2string
from pprint import pprint

class BondSetting(Setting):
    """
    BondSetting object

    The BondSetting object store the bondpair information.

    Parameters:

    label: str
        The label define the batoms object that a Bondsetting belong to.
    cutoff: float
        Cutoff used to calculate the maxmium bondlength for bond pairs.
    """
    
    def __init__(self, label, batoms = None, bondsetting = None, cutoff = 1.3) -> None:
        Setting.__init__(self, label, coll_name='%s_bond'%label)
        self.label = label
        self.name = 'bbond'
        self.cutoff = cutoff
        self.batoms = batoms
        if len(self) == 0:
            self.set_default(self.batoms.species.species_props, cutoff)
        if bondsetting is not None:
            for key, data in bondsetting.items():
                self[key] = data

    def build_materials(self, sp, select = 'sel0', node_inputs = None, 
                material_style = 'default'):
        """
        """
        from batoms.material import create_material
        colors = [sp.color1, sp.color2]
        for i in range(2):
            name = '%s_%s_%s'%(self.label, sp.name, i)
            if name not in bpy.data.materials:
                create_material(name,
                            color = colors[i],
                            node_inputs = node_inputs,
                            material_style = material_style,
                            backface_culling = True)
    
    def build_instancer(self, sp, select = 'sel0', vertices = 32, 
                        shape = 'CYLINDER', shade_smooth = True):
        # only build the needed one
        for order in [sp.order]:
            for style in [int(sp.style)]:
                name = 'bond_%s_%s_%s_%s'%(self.label, sp.name, order, style)
                radius = sp.width
                if name in bpy.data.objects:
                    obj = bpy.data.objects.get(name)
                    bpy.data.objects.remove(obj, do_unlink = True)
                self.cylinder(order = order, style = style, vertices = vertices, depth = 1, radius = 1)
                obj = bpy.context.view_layer.objects.active
                obj.name = name
                obj.data.name = name
                obj.batoms.batom.radius = radius
                #
                self.batoms.coll.children['%s_instancer'%self.label].objects.link(obj)
                # bpy.data.scenes['Scene'].objects.unlink(bb.obj)
                if shade_smooth:
                    bpy.ops.object.shade_smooth()
                obj.hide_set(True)
                #
                self.build_materials(sp, select)
                self.assign_materials(sp, order, style)
                # obj.data.materials.append(materials[data[0]])
                bpy.context.view_layer.update()
        return obj
    
    def cylinder(self, order = 1, style = 1, vertices = 32, depth = 1.0, radius = 1.0):
        """
        create cylinder and subdivde to 2 parts
        """
        radius = radius/order
        if style == 2:
            depth = depth/20
            radius = radius
        bpy.ops.mesh.primitive_cylinder_add(vertices = vertices, depth = depth, radius = radius)
        obj = bpy.context.view_layer.objects.active
        me = obj.data
        # select edges for subdivde
        n = len(me.vertices)
        vertices = np.zeros(n*3, dtype=np.float64)
        me.vertices.foreach_get('co', vertices)  
        vertices = vertices.reshape((n, 3))
        #
        me.update()
        m = len(me.edges)
        selects = np.zeros(m, dtype=bool)
        for i in range(m):
            v0 = me.edges[i].vertices[0]
            v1 = me.edges[i].vertices[1]
            center = (vertices[v0] + vertices[v1])/2
            if np.isclose(center[2], 0):
                selects[i] = True
        #
        bpy.ops.object.mode_set(mode = 'EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode = 'OBJECT')
        me.edges.foreach_set('select', selects)
        bpy.ops.object.mode_set(mode = 'EDIT')
        bpy.ops.mesh.subdivide(number_cuts = 1, smoothness = 0, fractal_along_normal = 0)
        bpy.ops.object.mode_set(mode = 'OBJECT')
        # order 
        mod = obj.modifiers.new('arrayOrder', 'ARRAY')
        mod.relative_offset_displace = (2, 0, 0)
        mod.count = order
        nv = len(obj.data.vertices)
        for i in range(nv):
            obj.data.vertices[i].co[0] -= 2*(order - 1)*radius
        bpy.ops.object.modifier_apply(modifier='Array')
        #
        if style == 2:
            mod = obj.modifiers.new('arrayStyle', 'ARRAY')
            mod.relative_offset_displace = (0, 0, 2)
            mod.count = 10
            nv = len(obj.data.vertices)
            for i in range(nv):
                obj.data.vertices[i].co[2] -= (10 - 1)*depth
        bpy.ops.object.modifier_apply(modifier='Array')
        return obj

    def assign_materials(self, sp, order, style = 0):
        # sort element by occu
        me = self.instancers[sp.name]['%s_%s'%(order, style)].data
        me.materials.clear()
        if style in [0, 2]:
            for i in range(2):
                me.materials.append(self.materials[sp.name][0])
        elif style == 1:
            for i in range(2):
                me.materials.append(self.materials[sp.name][i])
        # find the face index for sp1 and sp2
        #
        npoly = len(me.polygons)
        centers = np.zeros(npoly*3, dtype=np.float64)
        me.polygons.foreach_get('center', centers)  
        centers = centers.reshape((npoly, 3))
        material_indexs = np.zeros(npoly, dtype='int')
        index = np.where(centers[:, 2] < 0)[0]
        material_indexs[index] = 1
        me.polygons.foreach_set('material_index', material_indexs)

    @property
    def instancers(self):
        return self.get_instancers()
    
    def get_instancers(self):
        instancers = {}
        for sp in self:
            instancers[sp.name] = {}
            for order in [1, 2, 3]:
                for style in [0, 1, 2]:
                    name = 'bond_%s_%s_%s_%s'%(self.label, sp.name, order, style)
                    instancers[sp.name]['%s_%s'%(order, style)] = bpy.data.objects.get(name)
        return instancers
    
    @property
    def materials(self):
        return self.get_materials()
    
    def get_materials(self):
        materials = {}
        for sp in self:
            materials[sp.name] = {}
            for i in range(2):
                name = '%s_%s_%s'%(self.label, sp.name, i)
                mat = bpy.data.materials.get(name)
                materials[sp.name][i] = mat
        return materials

    def __setitem__(self, index, setdict):
        """
        Set properties
        """
        if isinstance(index, str):
            raise Exception("Bond index should be a tuple or list, \
    e.g. h2o.bondseeting[('O', 'H')]")
        name = tuple2string(index)
        subset = self.find(name)
        if subset is None:
            subset = self.collection.add()
        subset.name = name
        subset.species1 = index[0]
        subset.species2 = index[1]
        for key, value in setdict.items():
            setattr(subset, key, value)
        subset.label = self.label
        subset.flag = True
        self.build_instancer(self[name])

    def __getitem__(self, index):
        name = tuple2string(index)
        item = self.find(name)
        if item is None:
            raise Exception('%s not in %s setting'%(name, self.name))
        return item
        
    def set_default(self, species_props, cutoff = 1.3, self_interaction = True):
        """
        """
        for sel, data in species_props.items():
            bondtable = get_bondtable(self.label, sel, data, cutoff=cutoff, 
                        self_interaction = self_interaction)
            for key, value in bondtable.items():
                self[key] = value
    
    def extend(self, other):
        for b in other:
            self[(b.species1, b.species2)] = b.as_dict()
        # new
        species1 = set(self.batoms.species.species_props)
        species2 = set(other.batoms.species.species_props)
        nspecies1 = species1 - species2
        nspecies2 = species2 - species1
        for sp1 in nspecies1:
            for sp2 in nspecies2:
                species = {sp1: self.batoms.species.species_props[sp1], sp2: other.batoms.species.species_props[sp2]}
                self.set_default(species, self.cutoff, self_interaction = False)
    
    def add(self, bondpairs):
        if isinstance(bondpairs, tuple):
            bondpairs = [bondpairs]
        if isinstance(bondpairs[0], (str, int, float)):
            bondpairs = [bondpairs]
        for bondpair in bondpairs:
            species = {sp: self.batoms.species.species_props[sp] for sp in bondpair}
            maxlength = species[bondpair[0]]['radius'] + \
                        species[bondpair[1]]['radius']
            self[bondpair] = {'max': maxlength*self.cutoff, 
                            'color1':species[bondpair[0]]['color'],
                            'color2':species[bondpair[1]]['color'],
                            }
            self.set_default(species)
    
    def copy(self, label):
        object_mode()
        bondsetting = self.__class__(label)
        for b in self:
            bondsetting[(b.species1, b.species2)] = b.as_dict()
        return bondsetting
    
    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Bondpair    min     max   Search_bond    Polyhedra style\n'
        for b in self.collection:
            s += '{:10s}  {:4.3f}   {:4.3f}      {:10s}   {:10s}  {:4s} \n'.format(\
                b.name, b.min, b.max, str(b.search), str(b.polyhedra), b.style)
        s += '-'*60 + '\n'
        return s
    
    @property    
    def cutoff_dict(self):
        cutoff_dict = {}
        for b in self.collection:
            cutoff_dict[(b.species1, b.species2)] = [b.min, b.max]
        return cutoff_dict
    
    @property    
    def maxcutoff(self):
        maxcutoff = 2
        for bl in self.cutoff_dict.values():
            if bl[1] > maxcutoff:
                maxcutoff = bl[1]
        return maxcutoff
    
    def search_bond_list(self, atoms, bondlists0, offsets_skin0):
        bondlist1 = []
        offsets_skin1 = []
        bondlist2 = []
        offsets_skin2 = []
        speciesarray = np.array(atoms.arrays['species'])
        for b in self:
            if b.search == 1:
                temp = bondlists0[(speciesarray[bondlists0[:, 0]] == b.species1) 
                              & (speciesarray[bondlists0[:, 1]] == b.species2)]
                bondlist1.extend(temp)
                temp = offsets_skin0[(speciesarray[offsets_skin0[:, 0]] == b.species1)]
                offsets_skin1.extend(temp)
            elif b.search == 2:
                #  recursively if either sp1 or sp2
                temp = bondlists0[((speciesarray[bondlists0[:, 0]] == b.species1) 
                                  & (speciesarray[bondlists0[:, 1]] == b.species2))]
                bondlist2.extend(temp)
                temp1 = temp.copy()
                temp1[:, 0] = temp[:, 1]
                temp1[:, 1] = temp[:, 0]
                temp1[:, 2:] = -temp[:, 2:]
                bondlist2.extend(temp1)
                temp = offsets_skin0[(speciesarray[offsets_skin0[:, 0]] == b.species1)
                                   | (speciesarray[offsets_skin0[:, 0]] == b.species2)]
                offsets_skin2.extend(temp)
        return np.array(offsets_skin1), np.array(bondlist1), np.array(offsets_skin2), np.array(bondlist2)

    def draw_bonds(self, ):
        """Draw bonds.
        calculate bond in all farmes, and save
        get the max number of bonds for each pair
        draw the bonds
        add shape key
        the extra bonds use the find bond data.
        """
        
        from batoms.bond import Bbond
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
        bond_datas = {}
        tstart = time()
        for f in range(nframe):
            print('update bond: ', f)
            positions = frames[f]
            # if len(frames_boundary) > 0:
            #     positions_boundary = frames_boundary[f]
            #     positions = positions + positions_boundary
            # if len(frames_search) > 0:
            #     positions_search = frames_search[f]
            #     positions = positions + positions_search
            self.bondlist = build_bondlists(species, positions, 
                        self.batoms.cell, self.batoms.pbc, self.cutoff_dict)
            bond_kinds = calc_bond_data(species, positions, 
                        self.batoms.cell, size, self.bondlist, self,
                        arrays['model_style'])
            if f == 0:
                bond_datas = bond_kinds
            else:
                for kind, bond_data in bond_kinds.items():
                    if kind not in bond_datas:
                        bond_datas[kind] = bond_kinds[kind]
                    else:
                        bond_datas[kind]['positions'].extend(bond_data['positions'])
                        if len(bond_data['positions']) > 0:
                            if bond_datas[kind]['nposition'] < len(bond_data['positions'][0]):
                                bond_datas[kind]['nposition'] = len(bond_data['positions'][0])
        # print('calc bond: {0:10.2f} s'.format(time() - tstart))
            # nframe = len(bond_datas['positions'])
            # if nframe == 0: continue
            # change to same length
            # find max
            # nbs = [bond_data['positions'][i].shape[0] for i in range(nframe)]
            # nb_max = max(nbs)
            # frames_bond = np.zeros((nframe, nb_max, 7))
            # for i in range(nframe):
            #     frames_bond[i, 0:nbs[i], :] = bond_data['positions'][i]
            #     frames_bond[i, nbs[i]:, 0:3] = frames[i][0]
        if len(bond_datas) == 0:
            return
        bb = Bbond(self.batoms.label, species, bond_datas, 
                    # battr_inputs=bond_data['battr_inputs'],
                    batoms = self.batoms, 
                    bondsetting = self,
                    )
        # self.coll.objects.link(bb.obj)
        # bpy.data.collections['Collection'].objects.unlink(bb.obj)
        # bb.set_frames()
        # bpy.context.scene.frame_set(self.batoms.nframe)
        print('draw bond: {0:10.2f} s'.format(time() - tstart))

def get_bondtable(label, select, speciesdict, cutoff = 1.3, self_interaction = True):
    """
    """
    from batoms.data import default_bonds
    bondtable = {}
    for species1 in speciesdict:
        element1 = species1.split('_')[0]
        color1 = speciesdict[species1]['color']
        radius1 = cutoff * speciesdict[species1]['radius']
        for species2 in speciesdict:
            if not self_interaction and species1 == species2: continue
            element2 = species2.split('_')[0]
            element_pair = (element1, element2)
            if element_pair not in default_bonds: continue
            pair12 = (species1, species2)
            pair21 = (species2, species1)
            # only add bond once, except for hedrogen bond
            if pair21 in bondtable and ((element1 != 'H') and (element2 !='H')): continue
            color2 = speciesdict[species2]['color']
            radius2 = cutoff * speciesdict[species2]['radius']
            bondmax = radius1 + radius2
            bondtable[pair12] = {
                    'flag': True,
                    'label': label,
                    'select': select,
                    'species1': species1, 
                    'species2': species2, 
                    'min': 0.0, 
                    'max': bondmax,
                    'search': default_bonds[element_pair][0], 
                    'polyhedra': default_bonds[element_pair][1],
                    'color1': color1,
                    'color2': color2,
                    'width': 0.10,
                    'order': 1,
                    'order_offset': 0.15,
                    'style': '1',
                    }
    # special for hydrogen bond
    hbs = [('H', 'O'), ('H', 'N')]
    for pair in bondtable:
        key = (pair[0].split('_')[0], pair[1].split('_')[0])
        if key in hbs:
            bondtable[pair]['min'] = 1.2
            bondtable[pair]['max'] = 2.1
            bondtable[pair]['search'] = 0
            bondtable[pair]['color1'] = [0.1, 0.1, 0.1, 1.0]
            bondtable[pair]['width'] = 0.03
            bondtable[pair]['style'] = '2'
    return bondtable

def build_bondlists(species, positions, cell, pbc, cutoff):
    """
    The default bonds are stored in 'default_bonds'
    Get all pairs of bonding atoms
    remove_bonds
    """
    from batoms.neighborlist import primitive_neighbor_list, neighbor_kdtree
    if len(cutoff) == 0: return {}
    #
    tstart = time()
    # nli, nlj, nlS = primitive_neighbor_list('ijS', pbc,
    #                                cell,
    #                                positions, cutoff, species=species,
    #                                self_interaction=False,
    #                                max_nbins=1e6)
    nli, nlj, nlS = neighbor_kdtree('ijS', species, 
                positions, cell, pbc,
            cutoff, use_scaled_positions = False)
    # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
    bondlists = np.append(np.array([nli, nlj], dtype=int).T, np.array(nlS, dtype=int), axis = 1)
    bondlists1 = bondlists.copy()
    bondlists1[:, [0, 1]] = bondlists1[:, [1, 0]]
    bondlists2 = np.concatenate((bondlists, bondlists1), axis=0)
    np.unique(bondlists2)
    return bondlists

def calc_bond_data(speciesarray, positions, cell, radii, 
            bondlists, bondsetting,
            model_styles):
    """
    """
    from ase.data import chemical_symbols
    chemical_symbols = np.array(chemical_symbols)
    # properties
    nb =len(bondlists)
    orders = np.zeros(nb, dtype = int) # 1, 2, 3
    styles = np.zeros(nb, dtype = int) # 0, 1, 2, 3
    widths = np.ones(nb, dtype = float) # 0, 1, 2, 3
    species1 = np.ones(nb, dtype = 'U4') # 0, 1, 2, 3
    species2 = np.ones(nb, dtype = 'U4') # 0, 1, 2, 3
    model_styles = model_styles[bondlists[:, 0]]
    if nb == 0:
        datas = {
        'species1': species1,
        'species2': species2,
        'centers':np.zeros((0, 3)),
        'normals':np.zeros((0, 3)),
        'lengths':np.zeros((0, 3)),
        'widths': widths,
        'orders': orders,
        'styles': styles,
        'model_style':np.array([], dtype = int),
        }
        return datas
    tstart = time()
    #------------------------------------
    # offsets
    offsets = np.dot(bondlists[:, 2:5], cell)
    # bond vector and length
    vec = positions[bondlists[:, 0]] - (positions[bondlists[:, 1]] + offsets)
    length = np.linalg.norm(vec, axis = 1)
    nvec = vec/length[:, None]
    pos = [positions[bondlists[:, 0]] - nvec*radii[bondlists[:, 1]][:, None]*0.5,
            positions[bondlists[:, 1]] + offsets + nvec*radii[bondlists[:, 1]][:, None]*0.5]
    vec = pos[0] - pos[1]
    lengths = np.linalg.norm(vec, axis = 1) + 1e-8
    normals = vec/lengths[:, None]
    centers = (pos[0] + pos[1])/2.0
    #---------------------------------------------
    # model_styles = model_styles[bondlists[:, 0]]
    for b in bondsetting:
        spi = b.species1
        spj = b.species2
        indi = (speciesarray[bondlists[:, 0]] == spi)
        indj = (speciesarray[bondlists[:, 1]] == spj)
        ind = indi & indj
        orders[ind] = b.order
        styles[ind] = int(b.style)
        widths[ind] = b.width
        species1[ind] = spi
        species2[ind] = spj
    datas = {
        'species1':species1,
        'species2':species2,
        'centers':centers,
        'normals':normals,
        'lengths':lengths,
        'widths':widths,
        'orders':orders,
        'styles':styles,
        'model_styles':model_styles,
    }
    # print('datas: ', datas)
    return datas


if __name__ == "__main__":
    pass