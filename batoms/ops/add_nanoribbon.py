"""
Use ASE's build function
"""

import bpy
from bpy.types import Operator
from bpy.props import (StringProperty,
                       IntProperty,
                       FloatProperty,
                       BoolProperty,
                       )
from ase.build import graphene_nanoribbon
from batoms import Batoms


class BuildNanoribbon(Operator):
    bl_idname = "nano.nanoribbon_add"
    bl_label = "Add Nanoribbon"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add Nanoribbon")

    label: StringProperty(
        name="Label", default='graphene',
        description="Label")

    type: StringProperty(
        name="Type", default='zigzag',
        description=("The orientation of the ribbon. Must be either ‘zigzag’ or ‘armchair’."))

    n: IntProperty(
        name="n", default=3,
        min=1, soft_max=20,
        description=("The width of the nanoribbon. For armchair nanoribbons, this n may be half-integer to repeat by half a cell."))

    m: IntProperty(
        name="m", default=6,
        min=1, soft_max=20,
        description="The length of the nanoribbon.")

    C_H: FloatProperty(
        name="C_H", default=1.09,
        min=0, soft_max=15,
        description="Carbon-hydrogen bond length. Default: 1.09 Angstrom.")

    C_C: FloatProperty(
        name="C_C", default=1.42,
        min=0, soft_max=15,
        description="Carbon-carbon bond length. Default: 1.42 Angstrom.")

    vacuum: FloatProperty(
        name="Vacuum", default=5.0,
        min=0, soft_max=15,
        description="Amount of vacuum added to non-periodic directions, if present.")

    saturated: BoolProperty(
        name="Saturated", default=False,
        description="If true, hydrogen atoms are placed along the edge.")

    magnetic: BoolProperty(
        name="Magnetic", default=False,
        description="Make the edges magnetic.")

    initial_mag: FloatProperty(
        name="initial_mag", default=1.2,
        min=0, soft_max=15,
        description="Magnitude of magnetic moment if magnetic.")

    sheet: BoolProperty(
        name="Sheet", default=False,
        description=("If true, make an infinite sheet instead of a ribbon (default: False)"))

    main_element: StringProperty(
        name="Main_element", default='C',
        description="Main_element")

    saturate_element: StringProperty(
        name="Saturate_element", default='H',
        description="saturate_element")

    def execute(self, context):
        atoms = graphene_nanoribbon(
            n=self.n,
            m=self.m,
            type=self.type,
            saturated=self.saturated,
            C_H=self.C_H,
            C_C=self.C_C,
            vacuum=self.vacuum,
            magnetic=self.magnetic,
            initial_mag=self.initial_mag,
            sheet=self.sheet,
            main_element=self.main_element,
            saturate_element=self.saturate_element,
        )
        batoms = Batoms(label=self.label, from_ase=atoms)
        batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj
        self.report({"INFO"}, "Add nanoribbon {}".format(self.label))
        return {'FINISHED'}
