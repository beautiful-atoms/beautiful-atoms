"""
select expression

all
none


Logical
not Sele
sele1 and sele2


Properties
element O


Coordinates
z<5



Chemcial classses

solvent: water
hydrogens: 
metals


"""

import bpy
from batoms.base import Setting
from batoms.bondsetting import BondSetting, build_bondlists, calc_bond_data
from batoms.polyhedrasetting import PolyhedraSetting, build_polyhedralists
from batoms.isosurfacesetting import IsosurfaceSetting
from batoms.planesetting import PlaneSetting
from batoms.mssetting import MSsetting
from batoms.ribbon import Ribbon
import numpy as np
from batoms.butils import object_mode, show_index
from time import time


class Select():
    """

    sel = batoms.select[key]
    sel.model_style = 1

    """
    def __init__(self, label, name = 'sel0', batoms = None,
                bondsetting = None) -> None:
        self.batoms = batoms
        self.label = label
        self.name = name
        if bondsetting is None:
            self.bondsetting = BondSetting(self.label)

    @property
    def batom_dict(self):
        return self.get_batom_dict()
    
    def get_batom_dict(self):
        batom_dict = {}
        for obj in self.batoms.coll_atom.objects:
            if obj.batoms.batom.select == self.name:
                batom_dict[obj.name] = Batom(name = obj.name)
        return batom_dict

    @property
    def model_style(self):
        return self.get_model_style()
    
    @model_style.setter
    def model_style(self, model_style):
        self.set_model_style(model_style)
    
    def get_model_style(self):
        return int(self.batoms.selects.collection[self.name].model_style)
    
    def set_model_style(self, model_style):
        self.batoms.selects.collection[self.name].model_style = str(model_style)
        self.draw(draw_isosurface = False)
    
    @property
    def radius_style(self):
        return self.get_radius_style()
    
    @radius_style.setter
    def radius_style(self, radius_style):
        self.set_radius_style(radius_style)
    
    def get_radius_style(self):
        return self.batoms.selects.collection[self.name].radius_style
    
    def set_radius_style(self, radius_style):
        self.batoms.selects.collection[self.name].radius_style = radius_style
        self.draw(draw_isosurface = False)
    
    @property
    def show(self):
        return self.get_show()
    
    @show.setter
    def show(self, show):
        self.set_show(show)
    
    def get_show(self):
        return self.batoms.selects.collection[self.name].show
    
    def set_show(self, show):
        self.batoms.selects.collection[self.name].show = show
        self.draw(draw_isosurface = False)
    
    @property
    def mask(self):
        return self.get_mask()
    
    @mask.setter
    def mask(self, mask):
        self.set_mask(mask)
    
    def get_mask(self):
        mask = self.batoms.arrays[self.name]
        mask = np.array(mask, bool)
        return mask
    
    def set_mask(self, mask):
        self.batoms.set_arrays({self.name: mask})
    
    @property
    def scale(self):
        return self.get_scale()
    
    @scale.setter
    def scale(self, scale):
        self.set_scale(scale)
    
    def get_scale(self):
        scale = {}
        for coll in [self.batoms, self.batoms_boundary, self.batoms_search]:
            for batom in coll.values():
                scale[batom.species] = batom.scale
        return scale
    
    def set_scale(self, scale):
        """Set scale for all atoms

        scale: float
        """
        for coll in [self.batom_dict]:
            for batom in coll.values():
                batom.scale = scale
    
    def draw_bonds(self):
        """Draw bonds.
        calculate bond in all farmes, and save
        get the max number of bonds for each pair
        draw the bonds
        add shape key
        the extra bonds use the find bond data.
        """
        
        from batoms.butils import add_keyframe_visibility
        from batoms.bond import Bbond
        # if not self.bondlist:
        object_mode()
        self.batoms.clean_atoms_objects('bond')
        frames = self.batoms.get_frames(self.batoms.batoms)
        frames_boundary = self.batoms.get_frames(self.batoms.batoms_boundary)
        frames_search = self.batoms.get_frames(self.batoms.batoms_search)
        nframe = len(frames)
        bond_datas = {}
        tstart = time()
        mask = self.mask
        for f in range(nframe):
            print('update bond: ', f)
            atoms = frames[f][mask]
            if len(frames_boundary) > 0:
                atoms_boundary = frames_boundary[f][mask]
                atoms = atoms + atoms_boundary
            if len(frames_search) > 0:
                atoms_search = frames_search[f][mask]
                atoms = atoms + atoms_search
            atoms.pbc = False
            self.bondlist = build_bondlists(atoms, self.bondsetting.cutoff_dict)
            bond_kinds = calc_bond_data(atoms, self.batoms.species_props, self.bondlist, self.bondsetting)
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
        for species, bond_data in bond_datas.items():
            nframe = len(bond_data['positions'])
            if nframe == 0: continue
            # change to same length
            # find max
            nbs = [bond_data['positions'][i].shape[0] for i in range(nframe)]
            nb_max = max(nbs)
            frames_bond = np.zeros((nframe, nb_max, 7))
            for i in range(nframe):
                frames_bond[i, 0:nbs[i], :] = bond_data['positions'][i]
                frames_bond[i, nbs[i]:, 0:3] = frames[i].positions[0]
            bb = Bbond(self.label, species, frames_bond, 
                        color = bond_data['color'],
                        width = bond_data['width'],
                        segments = bond_data['segments'],
                        battr_inputs=bond_data['battr_inputs'])
            self.batoms.coll.children['%s_bond'%self.label].objects.link(bb.obj)
            bpy.data.collections['Collection'].objects.unlink(bb.obj)
            bb.set_frames()
        bpy.context.scene.frame_set(self.batoms.nframe)
        print('draw bond: {0:10.2f} s'.format(time() - tstart))

    def draw_polyhedras(self, bondlist = None):
        """Draw polyhedra.

        Parameters:

        bondlist: dict
        """
        object_mode()
        self.clean_atoms_objects('polyhedra')
        atoms = self.get_atoms_with_boundary(X = True)
        if bondlist is None:
            bondlist = build_bondlists(atoms, self.bondsetting.cutoff_dict)
        polyhedra_kinds = build_polyhedralists(atoms, bondlist, 
                          self.bondsetting, self.polyhedrasetting)
        for species, polyhedra_data in polyhedra_kinds.items():
            name = '%s_%s_%s'%(self.label, 'polyhedra', species)
            draw_surface_from_vertices(name, 
                                datas = polyhedra_data,
                                use_smooth = False,
                                coll = self.coll.children['%s_polyhedra'%self.label],
                                )
            if polyhedra_data['show_edge']:
                name = '%s_%s_%s'%(self.label, 'polyhedra_edge', species)
                draw_cylinder(name = name, 
                            datas = polyhedra_data['edge'],
                            coll = self.coll.children['%s_polyhedra'%self.label])
    
    def draw(self, model_style = None, draw_isosurface = True):
        """
        Draw atoms, bonds, polyhedra, .

        Parameters:

        model_style: str
        draw_isosurface: bool
        """
        if model_style is not None and model_style not in [0, 1, 2, 3]:
            raise Exception('model_style %s should be: 0, 1, 2, 3'%model_style)
        if not model_style:
            model_style = self.model_style
        else:
            self.model_style = model_style
        # self.draw_cell()
        bpy.ops.ed.undo_push()
        self.batoms.clean_atoms_objects('bond')
        bpy.ops.ed.undo_push()
        self.batoms.clean_atoms_objects('polyhedra')
        bpy.ops.ed.undo_push()
        if model_style == 0:
            self.scale = 1.0
        elif model_style == 1:
            self.scale = 0.4
            self.draw_bonds()
        elif model_style == 2:
            if self.polyhedra_style == 0:
                self.scale = 0.4
                self.draw_bonds()
                self.draw_polyhedras(self.bondlist)
            if self.polyhedra_style == 1:
                self.scale = 0.4
                self.draw_polyhedras()
            elif self.polyhedra_style == 2:
                self.scale = 0.01
                for b in self.bondsetting:
                    if b.polyhedra:
                        self.batoms[b.species1].scale = 0.4
                self.draw_polyhedras()
            elif self.polyhedra_style == 3:
                self.scale = 0.01
                self.draw_polyhedras()
        elif model_style == 3:
            self.scale = 0.01
            self.draw_bonds()

    def __repr__(self) -> str:
        s = '-'*60 + '\n'
        s = 'Name   model_style   radius_style     show \n'
        s += '{0:10s} {1:10s}   {2:10s}   {3:10s} \n'.format(\
                self.name, str(self.model_style), self.radius_style, str(self.show))
        s += '-'*60 + '\n'
        return s

class Selects(Setting):
    """
    sel = batoms.select.add(name, index)
    sel = batoms.select.from_dict({name: index})
    sel.model_style = 1

    """
    def __init__(self, label, batoms = None) -> None:
        self.label = label
        self.name = 'bselect'
        self.batoms = batoms

    @property
    def selects(self):
        return self.get_selects()
    
    def get_selects(self):
        selects = {}
        for b in self.collection:
            selects[b.name] = Select(self.label, name = b.name, batoms = self.batoms)
        return selects

    def __getitem__(self, key):
        """Return a subset of the Batom.

        key -- str, describing which batom to return.
        """
        if isinstance(key, str):
            if key not in self.selects:
                raise KeyError('%s is not in this structure'%key)
            return self.selects[key]
        elif isinstance(key, list):
            raise KeyError('dict not supported yet!')
    
    def __len__(self):
        n = 0
        for ba in self.selects.values():
            n += len(ba)
        return n
    
    def __repr__(self) -> str:
        s = '     Name        model_style    show  \n'
        i = 0
        for name, sel in self.selects.items():
            s += '{:3d}   {:10s}  {:5s}   {} \n'.format(i, 
                    name, str(sel.model_style), str(sel.show))
            i += 1
        return s
    
    def add(self, name, indices):
        data = {'name': name, 'indices': indices}
        self.from_dict(data)
        return self[name]

    def from_dict(self, data):
        """
        data = {'name': 'sel1', 'indices': [1, 2, 3], 'model_style': 1}
        """
        arrays = np.zeros(len(self.batoms), int)
        arrays[data['indices']] = 1
        arrays = {data['name']: arrays}
        data.pop('indices')
        self.batoms.set_arrays(arrays)
        self.batoms.split_select(data['name'])
        subset = self.find(data['name'])
        if subset is None:
            subset = self.collection.add()
        for key, value in data.items():
            setattr(subset, key, value)
        subset.label = self.label
        subset.flag = True
