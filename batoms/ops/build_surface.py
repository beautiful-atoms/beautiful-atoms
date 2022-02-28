"""
Use ASE's build function 
https://wiki.fysik.dtu.dk/ase/ase/build/surface.html?highlight=surfa#ase.build.surface
"""

from sympy import root
import bpy
from bpy.types import Operator
from bpy.props import (StringProperty,
                       IntProperty,
                       IntVectorProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       BoolProperty,
                       )
from ase.build import (fcc100, fcc110, fcc111, fcc111_root,
                       bcc100, bcc110, bcc111, bcc111_root,
                       hcp10m10, hcp0001, hcp0001_root,
                       diamond100, diamond111,
                       )
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
        min=0, soft_max=100,
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


class BuildSurfaceFCC111Root(Operator):
    bl_idname = "surface.fcc111_root"
    bl_label = "Add Root FCC(111) Surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add Root FCC(111) Surface")

    symbol: StringProperty(
        name="Symbol", default='Au',
        description="The chemical symbol of the element to use.")

    root: IntProperty(name="Root", min=1, soft_max=49, default=3)

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
        atoms = fcc111_root(self.symbol, size=self.size,
                       root=self.root,
                       a=a,
                       vacuum=self.vacuum,
                       orthogonal=self.orthogonal)
        Batoms(label=self.label, from_ase=atoms)
        return {'FINISHED'}


class BuildSurfaceBCC100(Operator):
    bl_idname = "surface.bcc100"
    bl_label = "Add BCC(100) Surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add BCC(100) Surface")

    symbol: StringProperty(
        name="Symbol", default='Fe',
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
        atoms = bcc100(self.symbol, size=self.size,
                       a=a,
                       vacuum=self.vacuum,
                       orthogonal=self.orthogonal,
                       periodic=self.periodic)
        Batoms(label=self.label, from_ase=atoms)
        return {'FINISHED'}


class BuildSurfaceBCC110(Operator):
    bl_idname = "surface.bcc110"
    bl_label = "Add BCC(110) Surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add BCC(110) Surface")

    symbol: StringProperty(
        name="Symbol", default='Fe',
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
        atoms = bcc110(self.symbol, size=self.size,
                       a=a,
                       vacuum=self.vacuum,
                       orthogonal=self.orthogonal,
                       periodic=self.periodic)
        Batoms(label=self.label, from_ase=atoms)
        return {'FINISHED'}


class BuildSurfaceBCC111(Operator):
    bl_idname = "surface.bcc111"
    bl_label = "Add BCC(111) Surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add BCC(111) Surface")

    symbol: StringProperty(
        name="Symbol", default='Fe',
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
        atoms = bcc111(self.symbol, size=self.size,
                       a=a,
                       vacuum=self.vacuum,
                       orthogonal=self.orthogonal,
                       periodic=self.periodic)
        Batoms(label=self.label, from_ase=atoms)
        return {'FINISHED'}


class BuildSurfaceBCC111Root(Operator):
    bl_idname = "surface.bcc111_root"
    bl_label = "Add Root BCC(111) Surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add Root BCC(111) Surface")

    symbol: StringProperty(
        name="Symbol", default='Fe',
        description="The chemical symbol of the element to use.")

    root: IntProperty(name="Root", min=1, soft_max=49, default=3)

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
        atoms = bcc111_root(self.symbol, size=self.size,
                       root=self.root,
                       a=a,
                       vacuum=self.vacuum,
                       orthogonal=self.orthogonal)
        Batoms(label=self.label, from_ase=atoms)
        return {'FINISHED'}

class BuildSurfaceHCP0001(Operator):
    bl_idname = "surface.hcp0001"
    bl_label = "Add HCP(0001) Surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add HCP(0001) Surface")

    symbol: StringProperty(
        name="Symbol", default='Ti',
        description="The chemical symbol of the element to use.")

    size: IntVectorProperty(
        name="Size", size=3, default=(1, 1, 4),
        min=1, soft_max=10,
        description="System size in units of the minimal unit cell.")

    a: FloatProperty(
        name="a", default=0,
        min=0, soft_max=100,
        description="Lattice constant.")

    c: FloatProperty(
        name="c", default=0,
        min=0, soft_max=100,
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
        c = None if self.c == 0 else self.c
        atoms = hcp0001(self.symbol, size=self.size,
                        a=a,
                        c=c,
                        vacuum=self.vacuum,
                        orthogonal=self.orthogonal,
                        periodic=self.periodic)
        Batoms(label=self.label, from_ase=atoms)
        return {'FINISHED'}


class BuildSurfaceHCP0001Root(Operator):
    bl_idname = "surface.hcp0001_root"
    bl_label = "Add Root HCP(0001) Surface"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add Root HCP(0001) Surface")

    symbol: StringProperty(
        name="Symbol", default='Au',
        description="The chemical symbol of the element to use.")

    root: IntProperty(name="Root", min=1, soft_max=49, default=3)

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
        atoms = hcp0001_root(self.symbol, size=self.size,
                       root=self.root,
                       a=a,
                       vacuum=self.vacuum,
                       orthogonal=self.orthogonal)
        Batoms(label=self.label, from_ase=atoms)
        return {'FINISHED'}
