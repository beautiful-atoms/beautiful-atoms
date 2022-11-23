"""
Use ASE's build function
https://wiki.fysik.dtu.dk/ase/ase/build/build.html?highlight=nanotube#ase.build.nanotube
"""

import bpy
from bpy.types import Operator
from bpy.props import (StringProperty,
                       IntProperty,
                       FloatProperty,
                       BoolProperty,
                       )
from ase.build import nanotube
from batoms import Batoms


class BuildNanotube(Operator):
    bl_idname = "nano.nanotube_add"
    bl_label = "Add Nanotube"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add Nanotube")

    label: StringProperty(
        name="Label", default='tube',
        description="Label")

    symbol: StringProperty(
        name="Symbol", default='C',
        description="The chemical symbol of the element to use.")

    n: IntProperty(
        name="n", default=6,
        min=1, soft_max=20,
        description="System size")

    m: IntProperty(
        name="m", default=0,
        min=0, soft_max=20,
        description="System size")

    length: IntProperty(
        name="Length", default=4,
        min=1, soft_max=20,
        description="System size")

    bond: FloatProperty(
        name="bond", default=1.4,
        min=0, soft_max=15,
        description="bond")

    vacuum: FloatProperty(
        name="vacuum", default=5.0,
        min=0, soft_max=15,
        description="vacuum")

    verbose: BoolProperty(
        name="verbose", default=True,
        description="verbose")

    def execute(self, context):
        if self.label == '':
            self.label = self.symbol
        atoms = nanotube(
            n=self.n,
            m=self.m,
            length=self.length,
            bond=self.bond,
            symbol=self.symbol,
            vacuum=self.vacuum,
            verbose=self.verbose)
        batoms = Batoms(label=self.label, from_ase=atoms)
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Add nanotube {}".format(self.label))
        return {'FINISHED'}
