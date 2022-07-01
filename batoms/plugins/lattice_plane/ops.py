import bpy
import bmesh
from bpy.props import (BoolProperty,
                       FloatProperty,
                       IntVectorProperty,
                       StringProperty
                       )
from batoms import Batoms
from batoms.ops.base import OperatorBatoms


class LatticePlaneAdd(OperatorBatoms):
    bl_idname = "plane.lattice_plane_add"
    bl_label = "Add Lattice Plane"
    bl_description = ("Add Lattice Plane to a Batoms")

    indices: IntVectorProperty(
        name="Miller indices", size=3, default=[0, 0, 1])

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=context.object.batoms.label)
        batoms.lattice_plane.settings.add(self.indices)
        context.view_layer.objects.active = obj
        self.report({"INFO"}, "Add lattice plane ({}{}{})."
            .format(self.indices[0], self.indices[1], self.indices[2]))
        return {'FINISHED'}


class LatticePlaneRemove(OperatorBatoms):
    bl_idname = "plane.lattice_plane_remove"
    bl_label = "Remove Lattice Plane"
    bl_description = ("Remove Lattice Plane to a Batoms")

    name: StringProperty(
        name="name", default='1-1-1',
        description="Name of Lattice Plane to be removed")

    all: BoolProperty(name="all",
                      default=False,
                      description="Remove all Lattice Planes")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.lattice_plane.settings.remove((self.name))
        context.view_layer.objects.active = obj
        return {'FINISHED'}


class LatticePlaneDraw(OperatorBatoms):
    bl_idname = "plane.lattice_plane_draw"
    bl_label = "Draw Lattice Plane"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = ("Draw Lattice Plane")

    name: StringProperty(
        name="name", default="ALL",
        description="Name of Lattice Plane to be drawed")

    def execute(self, context):
        obj = context.object
        batoms = Batoms(label=obj.batoms.label)
        batoms.lattice_plane.draw(self.name)
        context.view_layer.objects.active = batoms.obj
        return {'FINISHED'}


class LatticePlaneModify(OperatorBatoms):
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
        v = [s.index for s in bm.select_history if isinstance(
            s, bmesh.types.BMVert)]
        batoms = Batoms(label=obj.batoms.label)
        for i in v:
            setattr(batoms.bonds[i], self.key, getattr(self, self.key))
        # batoms.draw()
        context.view_layer.objects.active = obj
        return {'FINISHED'}
