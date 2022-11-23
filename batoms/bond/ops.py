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
        name="species", default='',
        description="Replaced by this species")

    species2: StringProperty(
        name="species", default='',
        description="Replaced by this species")

    def execute(self, context):
        if self.species1 == '' or self.species2 == '':
            return {'FINISHED'}
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        pair = (self.species1, self.species2)
        batoms.bond.settings.add(pair)
        context.view_layer.objects.active = obj
        self.report({"INFO"}, "Add bond pair {} {}".format(
            self.species1, self.species2))
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
        batoms.bond.settings.remove((self.name))
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
        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class BondOrderAutoSet(OperatorBatoms):
    bl_idname = "bond.bond_order_auto_set"
    bl_label = "Set bond order"
    bl_description = ("Set bond order to a Batoms")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        batoms.bond.bond_order_auto_set()
        context.view_layer.objects.active = obj
        return {'FINISHED'}

class BondShowHydrogenBond(OperatorBatoms):
    bl_idname = "bond.show_hydrogen_bond"
    bl_label = "Show hydrogen bond"
    bl_description = ("Show hydrogen bond.")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        batoms.bond.show_hydrogen_bond = not batoms.bond.show_hydrogen_bond
        context.view_layer.objects.active = obj
        return {'FINISHED'}

class BondShowSearch(OperatorBatoms):
    bl_idname = "bond.show_search"
    bl_label = "Show atoms by searching bonds"
    bl_description = ("Show atoms by searching bonds.")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        batoms.bond.show_search = not batoms.bond.show_search
        batoms.bond.update()
        context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}
