"""
"""
import collections
import bpy
from batoms.butils import object_mode
import numpy as np
from time import time
from pprint import pprint

def tuple2string(index):
    if isinstance(index, (str, int, float)):
        return str(index)
    name = str(index[0])
    for key in index[1:]:
        name += '-%s'%key
    return name

class Setting():
    """
    Setting object

    The Setting object store the other information for batoms.

    Parameters:

    label: str
        The label define the batoms object that a Setting belong to.
    """
    def __init__(self, label) -> None:
        self.label = label
        self.name = 'base'
    def get_data(self):
        data = {}
        for b in self.collection:
            data[b.name] = b
        return data
    @property
    def collection(self):
        return self.get_collection()
    def get_collection(self):
        collection = getattr(bpy.data.collections[self.label], self.name)
        return collection
    @property
    def species(self):
        return self.get_species()
    def get_species(self) -> dict:
        """
        read species from collection.
        """
        species = {}
        coll_atom = bpy.data.collections['%s_atom'%self.label]
        for ba in coll_atom.objects:
            species[ba.batom.species] = {
                'color': ba.children[0].data.materials[0].diffuse_color,
                'radius': ba.batom.radius}
        return species
    @property
    def data(self):
        return self.get_data()
    def __getitem__(self, index):
        name = tuple2string(index)
        item = self.find(name)
        if item is None:
            raise Exception('%s not in %s setting'%(name, self.name))
        return item
    def __setitem__(self, index, setdict):
        """
        Set properties
        """
        name = tuple2string(index)
        subset = self.find(name)
        if subset is None:
            subset = self.collection.add()
        subset.name = name
        for key, value in setdict.items():
            setattr(subset, key, value)
        subset.label = self.label
        subset.flag = True
    def copy(self, label):
        object_mode()
        setting = self.__class__(label)
        for key, b in self.data.items():
            setting[key] = b.as_dict()
        return setting
    def remove(self, index):
        """
        index: list of tuple
        """
        if isinstance(index, (str, int, float)):
            index = [(index)]
        if isinstance(index[0], (str, int, float)):
            index = [index]
        for key in index:
            name = tuple2string(key)
            i = self.collection.find(name)
            if i != -1:
                self.collection.remove(i)
            else:
                raise Exception('%s is not in %s'%(name, self.name))
    def __delitem__(self, index):
        self.remove(index)
    def __add__(self, other):
        self += other
        return self
    def __iadd__(self, other):
        self.extend(other)
        return self
    def extend(self, other):
        for key, value in other.data.items():
            self[key] = value.as_dict()
        return self
    def __iter__(self):
        item = self.collection
        for i in range(len(item)):
            yield item[i]
    def __len__(self):
        return len(self.collection)
    def find(self, name):
        i = self.collection.find(str(name))
        if i == -1:
            # print('%s is not in %s'%(name, self.name))
            return None
        else:
            return self.collection[i]


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
    def __init__(self, label, bondsetting = None, cutoff = 1.3) -> None:
        Setting.__init__(self, label)
        self.label = label
        self.name = 'bbond'
        self.cutoff = cutoff
        if len(self) == 0:
            self.set_default(self.species, cutoff)
        if bondsetting is not None:
            for key, data in bondsetting.items():
                self[key] = data
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
    def set_default(self, species, cutoff = 1.3, self_interaction = True):
        """
        """
        bondtable = get_bondtable(self.label, species, cutoff=cutoff, 
                    self_interaction = self_interaction)
        for key, value in bondtable.items():
            self[key] = value
    def extend(self, other):
        for b in other:
            self[(b.species1, b.species2)] = b.as_dict()
        # new
        species1 = set(self.species)
        species2 = set(other.species)
        nspecies1 = species1 - species2
        nspecies2 = species2 - species1
        for sp1 in nspecies1:
            for sp2 in nspecies2:
                species = {sp1: self.species[sp1], sp2: other.species[sp2]}
                self.set_default(species, self.cutoff, self_interaction = False)
    def add(self, bondpairs):
        if isinstance(bondpairs, tuple):
            bondpairs = [bondpairs]
        if isinstance(bondpairs[0], (str, int, float)):
            bondpairs = [bondpairs]
        for bondpair in bondpairs:
            species = {sp: self.species[sp] for sp in bondpair}
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
        s = 'Bondpair      min     max   Search_bond    Polyhedra style\n'
        for b in self.collection:
            s += '{0:10s} {1:4.3f}   {2:4.3f}      {3:10s}   {4:10s}  {5:4s} \n'.format(\
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


def get_bondtable(label, speciesdict, cutoff = 1.3, self_interaction = True):
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
                    'species1': species1, 
                    'species2': species2, 
                    'min': 0.5, 
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

def build_bondlists(atoms, cutoff):
    """
    The default bonds are stored in 'default_bonds'
    Get all pairs of bonding atoms
    remove_bonds
    """
    from batoms.neighborlist import neighbor_list
    if len(cutoff) == 0: return {}
    #
    tstart = time()
    nli, nlj, nlS = neighbor_list('ijS', atoms, cutoff, self_interaction=False)
    # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
    bondlists = np.append(np.array([nli, nlj], dtype=int).T, np.array(nlS, dtype=int), axis = 1)
    return bondlists

def calc_bond_data(atoms, species_props, bondlists, bondsetting):
    """
    """
    from ase.data import chemical_symbols
    positions = atoms.positions
    chemical_symbols = np.array(chemical_symbols)
    if 'species' not in atoms.arrays:
        atoms.new_array('species', np.array(atoms.get_chemical_symbols(), dtype = 'U20'))
    speciesarray = np.array(atoms.arrays['species'])
    bond_kinds = {}
    if len(bondlists) == 0:
        bond_kinds = {}
        return bond_kinds
    tstart = time()
    for b in bondsetting:
        spi = b.species1
        spj = b.species2
        indi = (speciesarray[bondlists[:, 0]] == spi)
        indj = (speciesarray[bondlists[:, 1]] == spj)
        bondlists1 = bondlists[indi & indj]
        if len(bondlists1) == 0: continue
        #------------------------------------
        # bond vector and length
        offset = bondlists1[:, 2:5]
        R = np.dot(offset, atoms.cell)
        vec = positions[bondlists1[:, 0]] - (positions[bondlists1[:, 1]] + R)
        length = np.linalg.norm(vec, axis = 1)
        nvec = vec/length[:, None]
        pos = [positions[bondlists1[:, 0]] - nvec*species_props[spi]['size']*0.5,
                positions[bondlists1[:, 1]] + R + nvec*species_props[spj]['size']*0.5]
        vec = pos[0] - pos[1]
        length = np.linalg.norm(vec, axis = 1) + 1e-8
        nvec = vec/length[:, None]
        center0 = (pos[0] + pos[1])/2.0
        # verts, faces, for instancing
        #---------------------------------------------
        # name of bond objects
        kinds = [('%s_%s'%(spi, spj), spi, b.color1[:]), 
                    ('%s_%s_%s'%(spi, spj, spi), spi, b.color1[:]), 
                    ('%s_%s_%s'%(spi, spj, spj), spj, b.color2[:]),
                ]
        for kind, species, color in kinds:
            if kind not in bond_kinds:
                battr = b.as_dict()
                battr.update({'species': species})
                bond_kinds[kind] = {'battr_inputs': {'bbond': battr}, 
                                    'species': species,
                                    'color': color,
                                    'width': b.width,
                                    'segments': b.segments,
                                    'style': b.style,
                                    'positions': [],
                                    'nposition': 0,
                                    'centers': [],
                                    'lengths': [],
                                    'normals': [],
                                    }
        #--------------------
        # bond order
        if b.order > 1:
            nbond = len(bondlists1)
            indi1 = (speciesarray[bondlists[:, 1]] == spi)
            indj2 = (speciesarray[bondlists[:, 0]] == spj)
            bondlists2 = bondlists[indi | indj | indi1 | indj2]
            high_order_offsets = []
            for i in range(nbond):
                # find another bond
                mask = np.logical_not(((bondlists2[:, 0] == bondlists1[i, 0]) 
                        & (bondlists2[:, 1] == bondlists1[i, 1]))
                        | ((bondlists2[:, 0] != bondlists1[i, 0])
                        & (bondlists2[:, 0] != bondlists1[i, 1]) 
                        & (bondlists2[:, 1] != bondlists1[i, 0])
                        & (bondlists2[:, 1] != bondlists1[i, 1])))
                localbondlist = bondlists2[mask]
                # determine the plane of high order bond
                if len(localbondlist) == 0:
                    second_bond = np.array([0.0, 0.0, 1])
                else:
                    second_bond = positions[localbondlist[0, 0]] - positions[localbondlist[0, 1]]
                norml = np.cross(second_bond, nvec[i]) + np.array([0.000000001, 0, 0])
                high_order_offset = np.cross(norml, nvec[i]) + np.array([0.000000001, 0, 0])
                high_order_offset = high_order_offset/np.linalg.norm(high_order_offset)
                high_order_offsets.append(high_order_offset)
            high_order_offsets = np.array(high_order_offsets)
            center1 = center0 + high_order_offsets*b.order_offset
            center2 = center0 - high_order_offsets*b.order_offset
            if b.order == 2:
                center0 = np.concatenate((center1, center2), axis = 0)
                nvec = np.tile(nvec, (2, 1))
                length = np.tile(length, 2)
            if b.order == 3:
                center0 = np.concatenate((center0, center1, center2), axis = 0)
                nvec = np.tile(nvec, (3, 1))
                length = np.tile(length, 3)
        # Unicolor cylinder
        if b.style == '0':
                bond_kinds[kinds[0][0]]['positions'] = [np.hstack((center0, nvec, length.reshape(-1, 1)))]
                bond_kinds[kinds[0][0]]['nposition'] = len(center0)
        # Bicolor cylinder
        elif b.style == '1':
            length = length/2.0
            for i in range(1, 3):
                center = center0 - (i - 1.5)*nvec*length[:, None]
                bond_kinds[kinds[i][0]]['positions'] = [np.hstack((center, nvec, length.reshape(-1, 1)))]
                bond_kinds[kinds[i][0]]['nposition'] = len(center)
        # Dashed line
        elif b.style == '2':
            length = length/2.0
            for i in range(1, 3):
                bond_kinds[kinds[i][0]]['vertices'] = 6
                step = 0.1
                maxlength = length.max()
                center = center0 - (i - 1.5)*nvec*step
                for d in np.arange(0, maxlength, step):
                    offset = 2*(i - 1.5)*nvec*d
                    center1 = center + offset
                    ind = np.where(d<length)[0]
                    bond_kinds[kinds[i][0]]['centers'].extend(center1[ind])
                    bond_kinds[kinds[i][0]]['lengths'].extend([step/2]*len(ind))
                    bond_kinds[kinds[i][0]]['normals'].extend(nvec[ind])
                bond_kinds[kinds[i][0]]['positions'] = [np.hstack((bond_kinds[kinds[i][0]]['centers'], 
                                bond_kinds[kinds[i][0]]['normals'], 
                                np.array(bond_kinds[kinds[i][0]]['lengths']).reshape(-1, 1)))]
                bond_kinds[kinds[i][0]]['nposition'] = len(bond_kinds[kinds[i][0]]['centers'])
        # Dotted line
        elif b.style == '3':
            length = length/2.0
            for i in range(1, 3):
                bond_kinds[kinds[i][0]]['vertices'] = 6
                step = 0.05
                maxlength = length.max()
                center = center0 - (i - 1.5)*nvec*step
                for d in np.arange(0, maxlength, step):
                    offset = 2*(i - 1.5)*nvec*d
                    center1 = center + offset
                    ind = np.where(d<length)[0]
                    bond_kinds[kinds[i][0]]['centers'].extend(center1[ind])
                    bond_kinds[kinds[i][0]]['lengths'].extend([step/4]*len(ind))
                    bond_kinds[kinds[i][0]]['normals'].extend(nvec[ind])
                bond_kinds[kinds[i][0]]['positions'] = [np.hstack((bond_kinds[kinds[i][0]]['centers'], 
                                bond_kinds[kinds[i][0]]['normals'], 
                                np.array(bond_kinds[kinds[i][0]]['lengths']).reshape(-1, 1)))]
                bond_kinds[kinds[i][0]]['nposition'] = len(bond_kinds[kinds[i][0]]['centers'])
    # print('calc_bond_data: {0:10.2f} s'.format(time() - tstart))
    return bond_kinds


if __name__ == "__main__":
    pass