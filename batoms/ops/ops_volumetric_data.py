import bpy
import bmesh
from bpy_extras.io_utils import ImportHelper, ExportHelper
from bpy.types import Operator
from bpy.props import (BoolProperty,
                       FloatProperty,
                       StringProperty,
                       EnumProperty
                       )
from batoms import Batoms
from batoms.ops.base import OperatorBatoms
from batoms.internal_data.bpy_data import get_volumetric_data



class VolumetricDataAdd(OperatorBatoms, ImportHelper):
    bl_idname = "batoms.volumetric_data_add"
    bl_label = "Add volumetric data"
    bl_description = ("Add volumetric data to a Batoms")

    filename_ext = ".cube"

    name: StringProperty(
        name="name", default='',
        description="name to be added")

    def execute(self, context):
        from ase.io.cube import read_cube_data
        import numpy as np
        if self.name == '':
            return {'FINISHED'}
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        if '.cube' in self.filepath:
            data, atoms = read_cube_data(bpy.path.abspath(self.filepath))
            cell = atoms.cell
        elif 'CHGCAR' in self.filepath:
            from pymatgen.io.vasp.outputs import VolumetricData
            poscar, data, data_aug = VolumetricData.parse_file(self.filepath)
            data = data['total']
            cell = poscar.structure.lattice.matrix
        # check demsion
        assert np.isclose(cell, batoms.cell).all()
        batoms.volumetric_data[self.name] = data
        context.view_layer.objects.active = obj
        self.report({"INFO"}, "Volumetric data {} is added.".format(self.name))
        return {'FINISHED'}


class VolumetricDataRemove(OperatorBatoms):
    bl_idname = "batoms.volumetric_data_remove"
    bl_label = "Remove volumetric data"
    bl_description = ("Remove volumetric data to a Batoms")

    name: StringProperty(
        name="name", default='C',
        description="name to be removed")

    all: BoolProperty(name="all",
                      default=False,
                      description="Remove all data")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.volumetric_data.remove((self.name))
        context.view_layer.objects.active = obj
        self.report({"INFO"}, "Volumetric data {} is removed.".format(self.name))
        return {'FINISHED'}

class VolumetricDataCreate(OperatorBatoms):
    bl_idname = "batoms.volumetric_data_create"
    bl_label = "Create volumetric data"
    bl_description = ("Create volumetric data to a Batoms. Add or minus.")

    name: StringProperty(
        name="name", default='',
        description="Volumetric data to be Created.")

    select_data1: EnumProperty(
        name="select_data1",
        items=get_volumetric_data,
        description="Volumetric data to be used.",
        default=0
    )

    select_data2: EnumProperty(
        name="select_data2",
        items=get_volumetric_data,
        description="Volumetric data to be used.",
        default=0
    )

    operator: EnumProperty(
        name="operator",
        items=[('Add', 'Add', '', 0),
                ('Minus', 'Minus', '', 1)
              ],
        description="Volumetric data to be used.",
        default=0,
        )


    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        if len(batoms.volumetric_data) < 2:
            self.report({"WARNING"}, "Not enough volumetric data to be used.".format(self.name))
            return {'FINISHED'}
        if self.name in batoms.volumetric_data.keys():
            self.name = "{}_1".format(self.name)
        if self.operator.upper() == 'ADD':
            batoms.volumetric_data[self.name] = batoms.volumetric_data[self.select_data1] + \
                    batoms.volumetric_data[self.select_data2]
        else:
            batoms.volumetric_data[self.name] = batoms.volumetric_data[self.select_data1] - \
                    batoms.volumetric_data[self.select_data2]
        context.view_layer.objects.active = obj
        self.report({"INFO"}, "Volumetric data {} is created.".format(self.name))
        return {'FINISHED'}
