"""
Use ASE's build function 
https://wiki.fysik.dtu.dk/ase/ase/build/build.html?highlight=nanotube#
"""
import bpy
from bpy.types import Operator
from bpy.props import (StringProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       BoolProperty,
                       )
from ase.build import molecule, bulk
from ase import Atoms
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
            remove_collection(self.label, keep_batom=False)
        else:
            batoms_list = read_batoms_list()
            for label in batoms_list:
                coll = bpy.data.collections.get(label)
                remove_collection(label, keep_batom=False)
        return {'FINISHED'}


class AddMolecule(Operator):
    bl_idname = "batoms.molecule_add"
    bl_label = "Add molecule"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add molecule")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    formula: StringProperty(
        name="Formula", default='CH4',
        description="formula")

    def execute(self, context):
        if self.label == '':
            self.label = self.formula
        atoms = molecule(self.formula)
        if self.label in bpy.data.collections:
            self.label = "%s.001"%self.label
        batoms = Batoms(label=self.label, from_ase=atoms)
        batoms.model_style = 1
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}


class AddBulk(Operator):
    bl_idname = "batoms.bulk_add"
    bl_label = "Add bulk"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add bulk")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    formula: StringProperty(
        name="Formula", default='Au',
        description="formula")

    crystalstructure: StringProperty(
        name="Crystalstructure", default='',
        description="Crystalstructure")

    latticeconstant: FloatVectorProperty(
        name="Lattice constant", default=(0, 0, 0),
        min=0, soft_max=100,
        description="Lattice constant")

    alpha: FloatProperty(
        name="Angle", default=0,
        min=0, soft_max=360,
        description="Angle in degrees for rhombohedral lattice.")

    covera: FloatProperty(
        name="Covera", default=1.6329931,
        min=0, soft_max=10,
        description="c/a ratio used for hcp. Default is ideal ratio: sqrt(8/3).")

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
                     a=a, b=b, c=c,
                     orthorhombic=self.orthorhombic,
                     cubic=self.cubic)
        if self.label in bpy.data.collections:
            self.label = "%s.001"%self.label
        batoms = Batoms(label=self.label, from_ase=atoms)
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}


class AddAtoms(Operator):
    bl_idname = "batoms.atoms_add"
    bl_label = "Add atoms"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add atoms")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    formula: StringProperty(
        name="Formula", default='C',
        description="formula")

    def execute(self, context):
        if self.label == '':
            self.label = self.formula
        atoms = Atoms(self.formula)
        if self.label in bpy.data.collections:
            self.label = "%s.001"%self.label
        batoms = Batoms(label=self.label, from_ase=atoms)
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}

