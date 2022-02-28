"""
Use ASE's build function 
https://wiki.fysik.dtu.dk/ase/ase/build/build.html?highlight=nanotube#
"""
import bpy
from bpy.types import Operator
from bpy.props import (StringProperty,
                       IntProperty,
                       IntVectorProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       BoolProperty,
                       )
from ase.build import molecule, bulk
from ase import Atoms
from batoms.utils.butils import get_selected_batoms
from batoms import Batoms


class deleteBatoms(Operator):
    bl_idname = "batoms.delete"
    bl_label = "Delete batoms"
    bl_description = ("Delete batoms")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    def execute(self, context):
        from batoms.utils.butils import remove_collection, read_batoms_list
        if self.label != '':
            coll = bpy.data.collections.get(self.label)
            remove_collection(self.label)
        else:
            batoms_list = read_batoms_list()
            for label in batoms_list:
                coll = bpy.data.collections.get(label)
                remove_collection(label)
        return {'FINISHED'}


class AddMolecule(Operator):
    bl_idname = "batoms.add_molecule"
    bl_label = "Add molecule"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add molecule")

    formula: StringProperty(
        name="Formula", default='H2O',
        description="formula")
    label: StringProperty(
        name="Label", default='',
        description="Label")

    def execute(self, context):
        if self.label == '':
            self.label = self.formula
        atoms = molecule(self.formula)
        Batoms(label=self.label, from_ase=atoms)
        return {'FINISHED'}


class AddBulk(Operator):
    bl_idname = "batoms.add_bulk"
    bl_label = "Add bulk"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add bulk")

    formula: StringProperty(
        name="Formula", default='Au',
        description="formula")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    crystalstructure: StringProperty(
        name="Crystalstructure", default='',
        description="Crystalstructure")

    latticeconstant: FloatVectorProperty(
        name="Lattice constant", default=(0, 0, 0),
        min = 0, soft_max = 100,
        description = "Lattice constant")
    
    alpha: FloatProperty(
        name="Angle", default=0,
        min = 0, soft_max = 360,
        description = "Angle in degrees for rhombohedral lattice.")
    
    covera: FloatProperty(
        name="Covera", default=1.6329931,
        min = 0, soft_max = 10,
        description = "c/a ratio used for hcp. Default is ideal ratio: sqrt(8/3).")
        
    orthorhombic: BoolProperty(
        name="Orthorhombic", default=False,
        description="Orthorhombic")

    cubic: BoolProperty(
        name="Cubic", default=False,
        description="Cubic")

    def execute(self, context):
        if self.label == '':
            self.label = self.formula
        a = None if self.latticeconstant[0] == 0 else \
            self.latticeconstant[0]
        b = None if self.latticeconstant[1] == 0 else \
            self.latticeconstant[1]
        c = None if self.latticeconstant[2] == 0 else \
            self.latticeconstant[2]
        crystalstructure = None if self.crystalstructure == '' \
            else self.crystalstructure
        atoms = bulk(self.formula, crystalstructure=crystalstructure,
                     a = a, b = b, c = c,
                     orthorhombic=self.orthorhombic,
                     cubic=self.cubic)
        Batoms(label=self.label, from_ase=atoms)
        return {'FINISHED'}


class AddAtoms(Operator):
    bl_idname = "batoms.add_atoms"
    bl_label = "Add atoms"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add atoms")

    formula: StringProperty(
        name="Formula", default='H',
        description="formula")
    label: StringProperty(
        name="Label", default='h',
        description="Label")

    def execute(self, context):
        atoms = Atoms(self.formula)
        Batoms(label=self.label, from_ase=atoms)
        return {'FINISHED'}


class AddSurface(Operator):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    bl_idname = "batoms.add_surface"
    bl_label = "Add surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add surface")

    label: StringProperty(
        name="Label", default='surf',
        description="Label")
    indices: IntVectorProperty(
        name="Miller indices", size=3, default=(1, 1, 1),
        soft_min=-10, soft_max=10,
        description="Miller indices for the plane")
    size: IntVectorProperty(
        name="Size", size=3, default=(1, 1, 1),
        min=1, soft_max=10,
        description="System size in units of the minimal unit cell.")
    vacuum: FloatProperty(
        name="vacuum", default=5.0,
        min=0, soft_max=15,
        description="vacuum")
    nlayer: IntProperty(name="layers", min=1, soft_max=10, default=4)
    termination: StringProperty(
        name="termination", default='',
        description="termination")

    def execute(self, context):
        from ase.build import surface
        from ase.build.surfaces_with_termination import surfaces_with_termination
        selected_batoms = self.selected_batoms
        if len(selected_batoms) != 1:
            raise Exception('Please select one structure')
        batoms = Batoms(selected_batoms[0])
        bulk = batoms.as_ase()
        if len(self.termination) == 0:
            atoms = surface(bulk, self.indices,
                            self.nlayer, self.vacuum)
        else:
            atoms = surfaces_with_termination(bulk, self.indices,
                                              self.nlayer, self.vacuum,
                                              termination=self.termination)
        batoms.hide = True
        label = batoms.label + ''.join(str(i) for i in self.indices)
        atoms = atoms*self.size
        batoms = Batoms(label=label, from_ase=atoms, movie=True)
        batoms.translate([2, 2, 2])
        return {'FINISHED'}



class AddRootSurface(Operator):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    bl_idname = "batoms.add_root_surface"
    bl_label = "Add root surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add root surface")

    label: StringProperty(
        name="Label", default='rootsurf',
        description="Label")
    root: IntProperty(name="root", min=1, soft_max=49, default=3)
    size: IntVectorProperty(
        name="Size", size=3, default=(1, 1, 1),
        min=1, soft_max=10,
        description="System size in units of the minimal unit cell.")

    def execute(self, context):
        from ase.build import root_surface
        selected_batoms = self.selected_batoms
        if len(selected_batoms) != 1:
            raise Exception('Please select one structure')
        batoms = Batoms(selected_batoms[0])
        surf = batoms.as_ase()
        atoms = root_surface(surf, self.root)
        batoms.hide = True
        label = batoms.label + '_root'
        atoms = atoms*self.size
        batoms = Batoms(label=label, from_ase=atoms, movie=True)
        batoms.translate([2, 2, 2])
        return {'FINISHED'}
