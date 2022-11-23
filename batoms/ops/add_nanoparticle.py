"""
Use ASE's build function
https://wiki.fysik.dtu.dk/ase/ase/cluster/cluster.html?highlight=cluster#module-ase.cluster
"""

import bpy
from bpy.types import Operator
from bpy.props import (StringProperty,
                       IntProperty,
                       BoolProperty,
                       )
from ase.cluster import Decahedron, Icosahedron, Octahedron, wulff_construction
from batoms import Batoms


class BuildDecahedron(Operator):
    bl_idname = "nano.decahedron_add"
    bl_label = "Add Decahedron"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add decahedron")

    symbol: StringProperty(
        name="Symbol", default='Au',
        description="The chemical symbol of the element to use.")

    p: IntProperty(
        name="p", default=5,
        min=1, soft_max=10,
        description="System size")

    q: IntProperty(
        name="q", default=2,
        min=1, soft_max=5,
        description="System size")

    r: IntProperty(
        name="r", default=0,
        min=0, soft_max=3,
        description="System size")

    latticeconstant: IntProperty(
        name="Lattice Constant", default=0,
        min=0, soft_max=20,
        description="Lattice Constant")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    def execute(self, context):
        if self.label == '':
            self.label = self.symbol
        latticeconstant = None if self.latticeconstant == 0 else \
            self.latticeconstant
        atoms = Decahedron(
            symbol=self.symbol,
            p=self.p,
            q=self.q,
            r=self.r,
            latticeconstant=latticeconstant)
        batoms = Batoms(label=self.label, from_ase=atoms)
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Add nanoparticle {}".format(self.label))
        return {'FINISHED'}


class BuildIcosahedron(Operator):
    bl_idname = "nano.icosahedron_add"
    bl_label = "Add Icosahedron"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add Icosahedron")

    symbol: StringProperty(
        name="Symbol", default='Au',
        description="The chemical symbol of the element to use.")

    noshells: IntProperty(
        name="noshells", default=5,
        min=1, soft_max=10,
        description="System size")

    latticeconstant: IntProperty(
        name="Lattice Constant", default=0,
        min=0, soft_max=20,
        description="Lattice Constant")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    def execute(self, context):
        if self.label == '':
            self.label = self.symbol
        latticeconstant = None if self.latticeconstant == 0 else \
            self.latticeconstant
        atoms = Icosahedron(
            symbol=self.symbol,
            noshells=self.noshells,
            latticeconstant=latticeconstant)
        batoms = Batoms(label=self.label, from_ase=atoms)
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Add nanoparticle {}".format(self.label))
        return {'FINISHED'}


class BuildOctahedron(Operator):
    bl_idname = "nano.octahedron_add"
    bl_label = "Add Octahedron"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add Octahedron")

    symbol: StringProperty(
        name="Symbol", default='Au',
        description="The chemical symbol of the element to use.")

    length: IntProperty(
        name="length", default=5,
        min=1, soft_max=10,
        description="Number of atoms on the square edges of the complete octahedron.")

    cutoff: IntProperty(
        name="cutoff", default=0,
        min=1, soft_max=10,
        description="Number of layers cut at each vertex.")

    latticeconstant: IntProperty(
        name="Lattice Constant", default=0,
        min=0, soft_max=20,
        description="Lattice Constant")

    alloy: BoolProperty(
        name="Alloy", default=False,
        description="If True the L1_2 structure is used.")

    label: StringProperty(
        name="Label", default='',
        description="Label")

    def execute(self, context):
        if self.label == '':
            self.label = self.symbol
        latticeconstant = None if self.latticeconstant == 0 else \
            self.latticeconstant
        atoms = Octahedron(
            symbol=self.symbol,
            cutoff=self.cutoff,
            length=self.length,
            latticeconstant=latticeconstant,
            alloy=self.alloy)
        batoms = Batoms(label=self.label, from_ase=atoms)
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Add nanoparticle {}".format(self.label))
        return {'FINISHED'}
