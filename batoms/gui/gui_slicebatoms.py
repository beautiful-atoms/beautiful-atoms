import bpy
from bpy.types import (
    Panel,
)
from bpy.props import (
    FloatProperty,
    BoolProperty,
)
# TODO: 4.2+ support
from .. import __package__ as batoms
from batoms.utils.butils import get_selected_vertices
from batoms.utils.attribute import get_mesh_attribute_bmesh, set_mesh_attribute_bmesh


class Batom_PT_prepare(Panel):
    bl_label = "Batom"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {"DEFAULT_CLOSED"}
    bl_category = "Batoms"
    bl_idname = "BATOMS_PT_Batom"

    # @classmethod
    # def poll(cls, context):
    #     obj = context.object
    #     if obj:
    #         return obj.batoms.type == 'BATOMS' and obj.mode != 'EDIT'
    #     else:
    #         return False

    def draw(self, context):
        layout = self.layout
        batom = context.scene.batoms.batom
        layout.prop(batom, "scale")
        layout.prop(batom, "show")
        layout.prop(batom, "bond")


# ---------------------------------------------------
def get_selected_batom():
    """Get the indices and the object of the selected vertices"""
    # bpy.ops.object.mode_set(mode="OBJECT")
    context = bpy.context
    if context.object and context.object.batoms.type != "OTHER":
        indices = get_selected_vertices(context.object)
        if len(indices) > 0:
            return context.object, indices
        else:
            return context.object, []
    bpy.ops.object.mode_set(mode="EDIT")
    return None, []


def set_object_attr(key, value):
    """set the attribute of the object"""
    # bpy.ops.object.mode_set(mode="OBJECT")
    obj, indices = get_selected_batom()
    if obj is not None and len(indices) > 0:
        for index in indices:
            set_mesh_attribute_bmesh(obj, key, value, index=index)
        # bpy.context.view_layer.objects.active = batom.obj
    bpy.ops.object.mode_set(mode="EDIT")


def get_scale(self):
    obj, indices = get_selected_batom()
    if obj is not None and len(indices) > 0:
        return get_mesh_attribute_bmesh(obj, "scale", index=indices[0])[0]
    else:
        return 0


def set_scale(self, value):
    self["scale"] = value
    set_object_attr("scale", value)


def get_show(self):
    obj, indices = get_selected_batom()
    if obj is not None and len(indices) > 0:
        return bool(get_mesh_attribute_bmesh(obj, "show", index=indices[0])[0])
    else:
        return 0


def set_show(self, value):
    self["show"] = value
    set_object_attr("show", value)


class BatomProperties(bpy.types.PropertyGroup):
    def Callback_modify_size(self, context):
        batom = bpy.context.scene.batoms.batom
        size = batom.size
        bpy.ops.batoms.batom_modify(key="size", size=size)

    def Callback_modify_bond(self, context):
        batom = bpy.context.scene.batoms.batom
        bond = batom.bond
        bpy.ops.batoms.batom_modify(key="bond", bond=bond)

    scale: FloatProperty(
        name="scale",
        default=0.6,
        min=0,
        soft_max=2,
        description="scale",
        get=get_scale,
        set=set_scale,
    )

    size: FloatProperty(
        name="size",
        default=1.5,
        min=0,
        soft_max=4,
        description="size",
        update=Callback_modify_size,
    )

    bond: BoolProperty(
        name="Bond", default=True, description="bond", update=Callback_modify_bond
    )

    show: BoolProperty(
        name="Show", default=True, description="show", get=get_show, set=set_show
    )
