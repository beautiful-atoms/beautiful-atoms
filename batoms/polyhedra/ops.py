import bpy
import bmesh
from bpy.types import Operator
from bpy.props import (BoolProperty,
                       FloatProperty,
                       StringProperty
                       )
from batoms import Batoms
from batoms.ops.base import OperatorBatoms


class PolyhedraAdd(OperatorBatoms):
    bl_idname = "batoms.polyhedra_add"
    bl_label = "Add polyhedra"
    bl_description = ("Add polyhedra to a Batoms")

    species: StringProperty(
        name="species", default='',
        description="species to be added")

    def execute(self, context):
        if self.species == '':
            return {'FINISHED'}
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        batoms.polyhedra.settings.add(self.species)
        batoms.coll.Bpolyhedra.ui_list_index = len(batoms.polyhedra.settings) - 1
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class PolyhedraRemove(OperatorBatoms):
    bl_idname = "batoms.polyhedra_remove"
    bl_label = "Remove polyhedra"
    bl_description = ("Remove polyhedra to a Batoms")

    species: StringProperty(
        name="species", default='C',
        description="species to be removed")

    all: BoolProperty(name="all",
                      default=False,
                      description="Remove all polyhedras")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        index = batoms.coll.Bpolyhedra.ui_list_index
        batoms.polyhedra.settings.remove((self.species))
        batoms.coll.Bpolyhedra.ui_list_index = min(max(0, index - 1),
                                                 len(batoms.polyhedra.settings) - 1)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class PolyhedraDraw(OperatorBatoms):
    bl_idname = "batoms.polyhedra_draw"
    bl_label = "Update polyhedra"
    bl_description = ("Update polyhedra to a Batoms")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.model_style = '2'
        batoms.draw()
        context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}


class PolyhedraModify(Operator):
    bl_idname = "batoms.polyhedra_modify"
    bl_label = "Modify polyhedra"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Modify polyhedra")

    key: StringProperty(
        name="key", default='style',
        description="Replaced by this polyhedra")

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
