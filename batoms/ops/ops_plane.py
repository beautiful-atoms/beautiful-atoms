import bpy
import bmesh
from bpy.types import (Panel,
                       Operator,
                       )
from bpy_extras.object_utils import AddObjectHelper
from bpy.props import (BoolProperty,
                       EnumProperty,
                       FloatProperty,
                       IntProperty,
                       IntVectorProperty,
                       StringProperty
                       )
from batoms import Batoms

class LatticePlaneAdd(Operator):
    bl_idname = "plane.lattice_plane_add"
    bl_label = "Add Lattice Plane"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add Lattice Plane to a Batoms")

    indices: IntVectorProperty(
        name="Miller indices", size=3, default=[0, 0, 1])

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
        batoms.lattice_plane.setting.add(self.indices)
        context.view_layer.objects.active = obj
        return {'FINISHED'}

class LatticePlaneRemove(Operator):
    bl_idname = "plane.lattice_plane_remove"
    bl_label = "Remove Lattice Plane"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Remove Lattice Plane to a Batoms")

    name: StringProperty(
        name="name", default='1-1-1',
        description="Name of Lattice Plane to be removed")
    
    all: BoolProperty(name="all",
                       default=False,
                       description="Remove all Lattice Planes")

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
        index = batoms.coll.batoms.latticeplane_index
        batoms.lattice_plane.setting.remove((self.name))
        batoms.coll.batoms.latticeplane_index = min(max(0, index - 1),
                len(batoms.lattice_plane.setting) - 1)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class LatticePlaneDraw(Operator):
    bl_idname = "plane.lattice_plane_draw"
    bl_label = "Draw Lattice Plane"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Draw Lattice Plane")

    name: StringProperty(
        name="name", default="ALL",
        description="Name of Lattice Plane to be drawed")

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
        batoms.lattice_plane.draw(self.name)
        context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}

class LatticePlaneModify(Operator):
    bl_idname = "plane.lattice_plane_modify"
    bl_label = "Modify lattice_plane"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Modify lattice_plane")

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



class CrystalShapeAdd(Operator):
    bl_idname = "plane.crystal_shape_add"
    bl_label = "Add Lattice Plane"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Add Lattice Plane to a Batoms")

    indices: IntVectorProperty(
        name="Miller indices", size=3, default=[0, 0, 1])

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
        batoms.crystal_shape.setting.add(self.indices)
        context.view_layer.objects.active = obj
        return {'FINISHED'}

class CrystalShapeRemove(Operator):
    bl_idname = "plane.crystal_shape_remove"
    bl_label = "Remove Lattice Plane"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Remove Lattice Plane to a Batoms")

    name: StringProperty(
        name="name", default='1-1-1',
        description="Name of Lattice Plane to be removed")
    
    all: BoolProperty(name="all",
                       default=False,
                       description="Remove all Lattice Planes")

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
        index = batoms.coll.batoms.crystalshape_index
        batoms.crystal_shape.setting.remove((self.name))
        batoms.coll.batoms.crystalshape_index = min(max(0, index - 1),
                len(batoms.crystal_shape.setting) - 1)
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class CrystalShapeDraw(Operator):
    bl_idname = "plane.crystal_shape_draw"
    bl_label = "Draw Crystal Shape"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Draw Crystal Shape")

    name: StringProperty(
        name="name", default="ALL",
        description="Name of Crystal Shape to be drawed")

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
        batoms.crystal_shape.draw(self.name)
        context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}

class CrystalShapeModify(Operator):
    bl_idname = "plane.crystal_shape_modify"
    bl_label = "Modify lattice_plane"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Modify lattice_plane")

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
