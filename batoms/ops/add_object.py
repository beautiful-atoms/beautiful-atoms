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
            self.label = "%s_001" % self.label
        batoms = Batoms(label=self.label, from_ase=atoms)
        batoms.model_style = 1
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Add molecule {}".format(self.label))
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
            self.label = "%s_001" % self.label
        batoms = Batoms(label=self.label, from_ase=atoms)
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Add bulk {}".format(self.label))
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
            self.label = "%s_001" % self.label
        batoms = Batoms(label=self.label, from_ase=atoms)
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Add atoms {}".format(self.label))
        return {'FINISHED'}


class AddSmiles(Operator):
    bl_idname = "batoms.smiles_add"
    bl_label = "Add Smiles"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add Smiles")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    smiles: StringProperty(
        name="smiles", default='CCO',
        description="smiles")

    def execute(self, context):
        from openbabel import pybel
        if self.label == '':
            self.label = self.smiles
        if self.label in bpy.data.collections:
            self.label = "%s_001" % self.label
        mol = pybel.readstring("smi", self.smiles)
        mol.make3D(forcefield='mmff94', steps=100)
        batoms = Batoms(label=self.label, from_pybel=mol)
        batoms.model_style = 1
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Add smiles {}".format(self.label))
        return {'FINISHED'}
