"""
"""
import collections
import bpy
from batoms.butils import object_mode
import numpy as np
from time import time
from pprint import pprint


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
        item = self.find(index)
        if item is None:
            raise Exception('%s not in %s setting'%(index, self.name))
        return item
    def __setitem__(self, index, value):
        """
        Add bondpair one by one
        """
        collection = self.collection
        item = self.find(index)
        
    def copy(self, label):
        object_mode()
        bondsetting = self.__class__(label)
        for key, b in self.data.items():
            bondsetting[key] = b.as_list()
        return bondsetting
    def __add__(self, other):
        self += other
        return self
    def __iadd__(self, other):
        self.extend(other)
        return self
    def extend(self, other):
        for key, value in other.data.items():
            self[key] = value.as_list()
        return self
    def __iter__(self):
        item = self.collection
        for i in range(len(item)):
            yield item[i]
    def __len__(self):
        return len(self.collection)
    def find(self, name):
        i = self.collection.find(name)
        if i == -1:
            return None
        else:
            return self.collection[i]


class Bondsetting(Setting):
    """
    Bondsetting object

    The Bondsetting object store the bondpair information.

    Parameters:

    label: str
        The label define the batoms object that a Bondsetting belong to.
    cutoff: float
        Cutoff used to calculate the maxmium bondlength for bond pairs.
    """
    def __init__(self, label, cutoff = 1.3) -> None:
        Setting.__init__(self, label)
        self.label = label
        self.name = 'bbond'
        self.cutoff = cutoff
        self.set_default(self.species, cutoff)
    def __setitem__(self, index, value):
        """
        Add bondpair one by one
        """
        bond = self.find(index)
        if bond is None:
            bond = self.collection.add()
        bond.symbol1 = index.split('-')[0]
        bond.symbol2 = index.split('-')[1]
        bond.name = index
        bond.min = value[0]
        bond.max = value[1]
        bond.search = value[2]
        bond.polyhedra = value[3]
        if len(value) == 8:
            bond.color1 = value[4]
            bond.color2 = value[5]
            bond.width = value[6]
            bond.style = value[7]
    def set_default(self, species, cutoff = 1.3):
        """
        """
        bondtable = get_bondtable(species, cutoff=cutoff)
        for key, value in bondtable.items():
            self['%s-%s'%(key[0], key[1])] = value
    def extend(self, other):
        for key, value in other.data.items():
            self[key] = value.as_list()
        # new
        species1 = set(self.species)
        species2 = set(other.species)
        nspecies1 = species1 - species2
        nspecies2 = species2 - species1
        for sp1 in nspecies1:
            for sp2 in nspecies2:
                species = {sp1: self.species[sp1], sp2: other.species[sp2]}
                self.set_default(species, self.cutoff)
    def add_bonds(self, bondpair):
        for key in bondpair:
            species = {sp: self.species[sp] for sp in key}
            self.set_default(species)
    def remove_bonds(self, bondpair):
        if isinstance(bondpair[0], str):
            bondpair = [bondpair]
        for key in bondpair:
            name = '%s-%s'%(key[0], key[1])
            i = self.collection.find(name)
            if i != -1:
                self.collection.remove(i)
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
            cutoff_dict[(b.symbol1, b.symbol2)] = [b.min, b.max]
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
        speciesarray = np.array(atoms.info['species'])
        for b in self:
            if b.search == 1:
                temp = bondlists0[(speciesarray[bondlists0[:, 0]] == b.symbol1) 
                              & (speciesarray[bondlists0[:, 1]] == b.symbol2)]
                bondlist1.extend(temp)
                temp = offsets_skin0[(speciesarray[offsets_skin0[:, 0]] == b.symbol1)]
                offsets_skin1.extend(temp)
            elif b.search == 2:
                #  recursively if either sp1 or sp2
                temp = bondlists0[((speciesarray[bondlists0[:, 0]] == b.symbol1) 
                                  & (speciesarray[bondlists0[:, 1]] == b.symbol2))]
                bondlist2.extend(temp)
                temp1 = temp.copy()
                temp1[:, 0] = temp[:, 1]
                temp1[:, 1] = temp[:, 0]
                temp1[:, 2:] = -temp[:, 2:]
                bondlist2.extend(temp1)
                temp = offsets_skin0[(speciesarray[offsets_skin0[:, 0]] == b.symbol1)
                                   | (speciesarray[offsets_skin0[:, 0]] == b.symbol2)]
                offsets_skin2.extend(temp)
        return np.array(offsets_skin1), np.array(bondlist1), np.array(offsets_skin2), np.array(bondlist2)


def get_bondtable(speciesdict, cutoff = 1.3):
    """
    """
    from batoms.data import default_bonds
    bondtable = {}
    for species1 in speciesdict:
        element1 = species1.split('_')[0]
        color1 = speciesdict[species1]['color']
        radius1 = cutoff * speciesdict[species1]['radius']
        for species2 in speciesdict:
            element2 = species2.split('_')[0]
            pair = (element1, element2)
            if pair not in default_bonds: continue
            color2 = speciesdict[species2]['color']
            radius2 = cutoff * speciesdict[species2]['radius']
            bondmax = radius1 + radius2
            bondtable[(species1, species2)] = [0.5, bondmax, default_bonds[pair][0], 
                            default_bonds[pair][1], color1, color2, 0.10, '1']
    # special for hydrogen bond
    if ('H', 'O') in bondtable:
        bondtable[('H', 'O')] = [1.2, 2.1, 0, False, [0.1, 0.1, 0.1, 1.0], [0.1, 0.1, 0.1, 1.0], 0.03, '2']
    if ('H', 'N') in bondtable:
        bondtable[('H', 'N')] = [1.2, 2.1, 0, False, [0.1, 0.1, 0.1, 1.0], [0.1, 0.1, 0.1, 1.0], 0.03, '2']
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
    bondlists = np.append(np.array([nli, nlj], dtype=int).T, np.array(nlS, dtype=int), axis = 1)
    # print('build_bondlists: {0:10.2f} s'.format(time() - tstart))
    return bondlists

def calc_bond_data(batoms, bondlists, bondsetting):
    """
    """
    from ase.data import chemical_symbols
    atoms = batoms.get_atoms_with_boundary()
    positions = atoms.positions
    chemical_symbols = np.array(chemical_symbols)
    if 'species' not in atoms.info:
        atoms.info['species'] = atoms.get_chemical_symbols()
    speciesarray = np.array(atoms.info['species'])
    bond_kinds = {}
    if len(bondlists) == 0:
        bond_kinds = {}
        return
    tstart = time()
    for b in bondsetting:
        spi = b.symbol1
        spj = b.symbol2
        bondlists1 = bondlists[(speciesarray[bondlists[:, 0]] == spi) 
                    & (speciesarray[bondlists[:, 1]] == spj)]
        if len(bondlists1) == 0: continue
        offset = bondlists1[:, 2:5]
        R = np.dot(offset, atoms.cell)
        vec = positions[bondlists1[:, 0]] - (positions[bondlists1[:, 1]] + R)
        length = np.linalg.norm(vec, axis = 1)
        nvec = vec/length[:, None]
        pos = [positions[bondlists1[:, 0]] - nvec*batoms[spi].size*0.5,
                positions[bondlists1[:, 1]] + R + nvec*batoms[spj].size*0.5]
        vec = pos[0] - pos[1]
        length = np.linalg.norm(vec, axis = 1)
        nvec = vec/length[:, None]
        nvec = nvec + 1e-8
        # verts, faces, for instancing
        # v1 = nvec + np.array([1.2323, 0.493749, 0.5604937284])
        # tempv = np.einsum("ij, ij->i", v1, nvec)
        # v11 = v1 - (nvec.T*tempv).T
        # templengh = np.linalg.norm(v11, axis = 1)
        # v11 = v11/templengh[:, None]/2.828427
        # tempv = np.cross(nvec, v11)
        # v22 = (tempv.T*(length*length)).T
        #
        kinds = [('%s_%s'%(spi, spj), b.color1), 
                    ('%s_%s_%s'%(spi, spj, spi), b.color1), 
                    ('%s_%s_%s'%(spi, spj, spj), b.color2),
                ]
        for kind, color in kinds:
            if kind not in bond_kinds:
                bond_kinds[kind] = {'color': color[:3], 'verts': [], 'transmit': color[3], 
                                    'width': b.width,
                                    'centers': [],
                                    'lengths': [],
                                    'normals': [], 
                                    'style': b.style}
        center0 = (pos[0] + pos[1])/2.0
        # Unicolor cylinder
        if b.style == '0':
                bond_kinds[kinds[0][0]]['centers'] = center0
                bond_kinds[kinds[0][0]]['lengths'] = length
                bond_kinds[kinds[0][0]]['normals'] = nvec
        # Bicolor cylinder
        elif b.style == '1':
            length = length/2.0
            for i in range(1, 3):
                center = (center0 + pos[i - 1])/2.0
                bond_kinds[kinds[i][0]]['centers'] = center
                bond_kinds[kinds[i][0]]['lengths'] = length
                bond_kinds[kinds[i][0]]['normals'] = nvec
                # bond_kinds[kinds[i][0]]['verts'] = center + v11
                # bond_kinds[kinds[i][0]]['verts'] = np.append(bond_kinds[kinds[i][0]]['verts'], center - v11, axis = 0)
                # bond_kinds[kinds[i][0]]['verts'] = np.append(bond_kinds[kinds[i][0]]['verts'], center + v22, axis = 0)
                # bond_kinds[kinds[i][0]]['verts'] = np.append(bond_kinds[kinds[i][0]]['verts'], center - v22, axis = 0)
        # Dashed line
        elif b.style == '2':
            length = length/4.0
            for i in range(1, 3):
                bond_kinds[kinds[i][0]]['centers'] = []
                bond_kinds[kinds[i][0]]['lengths'] = []
                bond_kinds[kinds[i][0]]['normals'] = []
                step = 0.2
                maxlength = length.max()
                center = (center0 + pos[i - 1])/2.0
                for d in np.arange(-maxlength, maxlength, step):
                    # offset0 = np.linspace(-2, 2, nc)
                    # np.where(offset0>-length/2.0 & offset<length/2.0)
                    offset = nvec*d
                    center1 = center + offset
                    ind = np.where( (d>-length) & (d<length))[0]
                    bond_kinds[kinds[i][0]]['centers'].extend(center1[ind])
                    bond_kinds[kinds[i][0]]['lengths'].extend([step/2]*len(ind))
                    bond_kinds[kinds[i][0]]['normals'].extend(nvec[ind])
        # Dotted line
        elif b.style == '3':
            length = length/4.0
            for i in range(1, 3):
                bond_kinds[kinds[i][0]]['centers'] = []
                bond_kinds[kinds[i][0]]['lengths'] = []
                bond_kinds[kinds[i][0]]['normals'] = []
                step = 0.05
                maxlength = length.max()
                center = (center0 + pos[i - 1])/2.0
                for d in np.arange(-maxlength, maxlength, step):
                    # offset0 = np.linspace(-2, 2, nc)
                    # np.where(offset0>-length/2.0 & offset<length/2.0)
                    offset = nvec*d
                    center1 = center + offset
                    ind = np.where( (d>-length) & (d<length))[0]
                    bond_kinds[kinds[i][0]]['centers'].extend(center1[ind])
                    bond_kinds[kinds[i][0]]['lengths'].extend([step/4]*len(ind))
                    bond_kinds[kinds[i][0]]['normals'].extend(nvec[ind])
    # print('calc_bond_data: {0:10.2f} s'.format(time() - tstart))
    return bond_kinds





if __name__ == "__main__":
    from ase.build import bulk, molecule
    from ase.atoms import Atoms
    from ase.io import read
    from ase.visualize import view
    from ase.neighborlist import neighbor_list
    from boundary import search_boundary
    # atoms = bulk('Pt', cubic = True)
    atoms = read('docs/source/_static/datas/tio2.cif')
    # atoms = read('docs/source/_static/datas/mof-5.cif')
    atoms = atoms*[4, 4, 4]
    # positions1, offsets1, positions2, offsets2 = search_boundary(atoms.positions, atoms.cell, boundary=[[-0.6, 1.6], [-0.6, 1.6], [-0.6, 1.6]])
    bondsetting = {('Ti', 'O'): [0, 2.5, True, False], ('O', 'O'): [0, 1.5, False, False]}
    bondlists = build_bondlists(atoms, bondsetting)
    