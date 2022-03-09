import bpy
import bmesh
from bpy.types import (Panel,
                       Operator,
                       )
from bpy_extras.object_utils import AddObjectHelper
from bpy.props import (BoolProperty,
                       FloatProperty,
                       IntProperty,
                       StringProperty,
                       )
from batoms import Batoms

class BondPairAdd(Operator):
    bl_idname = "bond.bond_pair_add"
    bl_label = "Add Bond Pair"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add Bond Pair to a Batoms")

    species1: StringProperty(
        name="species", default='C',
        description="Replaced by this species")

    species2: StringProperty(
        name="species", default='H',
        description="Replaced by this species")

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            return obj.batoms.type != 'OTHER'
        else:
            return False

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        pair = (self.species1, self.species2)
        batoms.bonds.setting.add(pair)
        context.view_layer.objects.active = obj
        return {'FINISHED'}

class BondPairRemove(Operator):
    bl_idname = "bond.bond_pair_remove"
    bl_label = "Remove Bond Pair"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Remove Bond Pair to a Batoms")

    name: StringProperty(
        name="name", default='C-H',
        description="Name of bond pair to be removed")
    
    all: BoolProperty(name="all",
                       default=False,
                       description="Remove all bond pairs")

    @classmethod
    def poll(cls, context):
        obj = context.object
        if obj:
            return obj.batoms.type != 'OTHER'
        else:
            return False

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        index = batoms.coll.batoms.bond_index
        batoms.bonds.setting.remove((self.name))
        batoms.coll.batoms.bond_index = min(max(0, index - 1),
                len(batoms.bonds.setting) - 1)
        context.view_layer.objects.active = obj
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

    search: IntProperty(name="Search mode", default=0, 
                )
    polyhedra: BoolProperty(name="polyhedra", default=False, 
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
        v = [s.index for s in bm.select_history if isinstance(s, bmesh.types.BMVert)]
        batoms = Batoms(label=obj.batoms.label)
        for i in v:
            setattr(batoms.bonds[i], self.key, getattr(self, self.key))
        # batoms.draw()
        context.view_layer.objects.active = obj
        return {'FINISHED'}
