"""
Use ASE's build function 
https://wiki.fysik.dtu.dk/ase/ase/build/surface.html?highlight=surfa#ase.build.surface
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
from ase.build import molecule, bulk, fcc100, fcc110, fcc111
from ase import Atoms
from batoms.utils.butils import get_selected_batoms
from batoms import Batoms


class BuildSurfaceFCC100(Operator):
    bl_idname = "surface.fcc100"
    bl_label = "Add FCC(100) Surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add FCC(100) Surface")

    symbol: StringProperty(
        name="Symbol", default='Au',
        description="The chemical symbol of the element to use.")

    size: IntVectorProperty(
        name="Size", size=3, default=(1, 1, 4),
        min=1, soft_max=10,
        description="System size in units of the minimal unit cell.")

    a: FloatProperty(
        name="a", default=0,
        min=0, soft_max=100,
        description="Lattice constant.")

    vacuum: FloatProperty(
        name="vacuum", default=5.0,
        min=0, soft_max=15,
        description="vacuum")

    orthogonal: BoolProperty(
        name="orthogonal", default=True,
        description="orthogonal")

    periodic: BoolProperty(
        name="Periodic", default=False,
        description="Periodic")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    def execute(self, context):
        if self.label == '':
            self.label = self.symbol
        a = None if self.a == 0 else self.a
        atoms = fcc100(self.symbol, size=self.size,
                       a=a,
                       vacuum=self.vacuum,
                       orthogonal=self.orthogonal,
                       periodic=self.periodic)
        Batoms(label=self.label, from_ase=atoms)
        return {'FINISHED'}


class BuildSurfaceFCC110(Operator):
    bl_idname = "surface.fcc110"
    bl_label = "Add FCC(110) Surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add FCC(110) Surface")

    symbol: StringProperty(
        name="Symbol", default='Au',
        description="The chemical symbol of the element to use.")

    size: IntVectorProperty(
        name="Size", size=3, default=(1, 1, 4),
        min=1, soft_max=10,
        description="System size in units of the minimal unit cell.")

    a: FloatProperty(
        name="a", default=0,
        min=0, soft_max=110,
        description="Lattice constant.")

    vacuum: FloatProperty(
        name="vacuum", default=5.0,
        min=0, soft_max=15,
        description="vacuum")

    orthogonal: BoolProperty(
        name="orthogonal", default=True,
        description="orthogonal")

    periodic: BoolProperty(
        name="Periodic", default=False,
        description="Periodic")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    def execute(self, context):
        if self.label == '':
            self.label = self.symbol
        a = None if self.a == 0 else self.a
        atoms = fcc110(self.symbol, size=self.size,
                       a=a,
                       vacuum=self.vacuum,
                       orthogonal=self.orthogonal,
                       periodic=self.periodic)
        Batoms(label=self.label, from_ase=atoms)
        return {'FINISHED'}


class BuildSurfaceFCC111(Operator):
    bl_idname = "surface.fcc111"
    bl_label = "Add FCC(111) Surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add FCC(111) Surface")

    symbol: StringProperty(
        name="Symbol", default='Au',
        description="The chemical symbol of the element to use.")

    size: IntVectorProperty(
        name="Size", size=3, default=(1, 1, 4),
        min=1, soft_max=10,
        description="System size in units of the minimal unit cell.")

    a: FloatProperty(
        name="a", default=0,
        min=0, soft_max=111,
        description="Lattice constant.")

    vacuum: FloatProperty(
        name="vacuum", default=5.0,
        min=0, soft_max=15,
        description="vacuum")

    orthogonal: BoolProperty(
        name="orthogonal", default=False,
        description="orthogonal")

    periodic: BoolProperty(
        name="Periodic", default=False,
        description="Periodic")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    def execute(self, context):
        if self.label == '':
            self.label = self.symbol
        a = None if self.a == 0 else self.a
        atoms = fcc111(self.symbol, size=self.size,
                       a=a,
                       vacuum=self.vacuum,
                       orthogonal=self.orthogonal,
                       periodic=self.periodic)
        Batoms(label=self.label, from_ase=atoms)
        return {'FINISHED'}
