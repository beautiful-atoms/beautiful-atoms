import bpy
import bmesh
from bpy.types import Operator
from bpy.props import (BoolProperty,
                       EnumProperty,
                       IntProperty,
                       StringProperty,
                       )
from batoms import Batoms
from batoms.ops.base import OperatorBatoms


class BondPairAdd(OperatorBatoms):
    bl_idname = "bond.bond_pair_add"
    bl_label = "Add Bond Pair"
    bl_description = ("Add Bond Pair to a Batoms")

    species1: StringProperty(
        name="species", default='C',
        description="Replaced by this species")

    species2: StringProperty(
        name="species", default='H',
        description="Replaced by this species")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        pair = (self.species1, self.species2)
        batoms.bonds.setting.add(pair)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class BondPairRemove(OperatorBatoms):
    bl_idname = "bond.bond_pair_remove"
    bl_label = "Remove Bond Pair"
    bl_description = ("Remove Bond Pair to a Batoms")

    name: StringProperty(
        name="name", default='C-H',
        description="Name of bond pair to be removed")

    all: BoolProperty(name="all",
                      default=False,
                      description="Remove all bond pairs")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        index = batoms.coll.batoms.bond_index
        batoms.bonds.setting.remove((self.name))
        batoms.coll.batoms.bond_index = min(max(0, index - 1),
                                            len(batoms.bonds.setting) - 1)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class BondDraw(OperatorBatoms):
    bl_idname = "bond.draw"
    bl_label = "Draw Bond"
    bl_description = ("Draw Bond")

    model_style: EnumProperty(
        name="model_style",
        description="Structural models",
        items=(('0', "Space-filling", "Use ball"),
               ('1', "Ball-and-stick", "Use ball and stick"),
               ('2', "Polyhedral", "Use polyhedral"),
               ('3', "Wireframe", "Use wireframe")),
        default='1')

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.draw(self.model_style)
        context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}


class BondModify(Operator):
    bl_idname = "bond.bond_modify"
    bl_label = "Modify bond"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Modify bond")

    key: StringProperty(
        name="key", default='style',
        description="Replaced by this species")

    style: StringProperty(
        name="style", default='0',
        description="style")

    show: BoolProperty(name="show", default=False,
                            )
    order: IntProperty(name="Bond order", default=1,
                       )

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            return obj.batoms.type == 'BOND' and obj.mode == 'EDIT'
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
        context.view_layer.objects.active = obj
        bpy.ops.object.mode_set(mode = 'EDIT')
        return {'FINISHED'}
