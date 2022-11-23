import bpy
import bmesh
from bpy.types import Operator
from bpy.props import (BoolProperty,
                       FloatProperty,
                       StringProperty
                       )
from batoms import Batoms
from batoms.ops.base import OperatorBatoms


class SpeciesAdd(OperatorBatoms):
    bl_idname = "batoms.species_add"
    bl_label = "Add Species"
    bl_description = ("Add Species to a Batoms")

    species: StringProperty(
        name="species", default='',
        description="Species to be added")

    def execute(self, context):
        if self.species == '':
            return {'FINISHED'}
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        batoms.species.add(self.species)
        context.view_layer.objects.active = obj
        self.report({"INFO"}, "{} is added.".format(self.species))
        return {'FINISHED'}


class SpeciesRemove(OperatorBatoms):
    bl_idname = "batoms.species_remove"
    bl_label = "Remove Species"
    bl_description = ("Remove Species to a Batoms")

    species: StringProperty(
        name="species", default='C',
        description="Species to be removed")

    all: BoolProperty(name="all",
                      default=False,
                      description="Remove all Speciess")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.species.remove((self.species))
        context.view_layer.objects.active = obj
        self.report({"INFO"}, "{} is removed.".format(self.species))
        return {'FINISHED'}


class SpeciesUpdate(OperatorBatoms):
    bl_idname = "batoms.species_update"
    bl_label = "Update Species"
    bl_description = ("Update Species to a Batoms")


    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.species.update()
        context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}


class SpeciesModify(Operator):
    bl_idname = "batoms.species_modify"
    bl_label = "Modify Species"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Modify Species")

    key: StringProperty(
        name="key", default='style',
        description="Replaced by this species")

    slice: BoolProperty(name="slice", default=False,
                        )
    boundary: BoolProperty(name="boundary", default=False,
                           )
    distance: FloatProperty(name="distance",
                            description="Distance from origin",
                            default=1)

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            return obj.batoms.type == 'MS' and obj.mode == 'EDIT'
        else:
            return False

    def execute(self, context):
        obj = context.object
        data = obj.data
        bm = bmesh.from_edit_mesh(data)
        v = [s.index for s in bm.select_history if isinstance(
            s, bmesh.types.BMVert)]
        batoms = Batoms(label=obj.batoms.label)
        for i in v:
            setattr(batoms.bonds[i], self.key, getattr(self, self.key))
        # batoms.draw()
        return {'FINISHED'}
