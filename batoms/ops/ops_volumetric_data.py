import bpy
import bmesh
from bpy_extras.io_utils import ImportHelper, ExportHelper
from bpy.types import Operator
from bpy.props import (BoolProperty,
                       FloatProperty,
                       StringProperty
                       )
from batoms import Batoms
from batoms.ops.base import OperatorBatoms


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
        data, atoms = read_cube_data(bpy.path.abspath(self.filepath))
        # check demsion
        assert np.isclose(atoms.cell, batoms.cell).all()
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



