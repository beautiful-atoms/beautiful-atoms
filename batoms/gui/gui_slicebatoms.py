import bpy
from bpy.types import (Panel,
                       )
from bpy.props import (FloatProperty,
                       FloatVectorProperty,
                       BoolProperty,
                       EnumProperty,
                       )
from batoms.utils.butils import get_selected_vertices
from batoms.slicebatoms import SliceBatoms
from batoms.batoms import Batoms

class Batom_PT_prepare(Panel):
    bl_label = "Batom"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
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
def get_active_batom():
    context = bpy.context
    if context.object and context.object.batoms.type != 'OTHER':
        indices = get_selected_vertices(context.object)
        if len(indices) > 0:
            batoms = Batoms(label=context.object.batoms.label)
            slicebatoms = batoms[indices]
            return slicebatoms
        else:
            return None
    return None

def set_batom_attr(key, value):
    """
    """
    batom = get_active_batom()
    if batom is not None:
        setattr(batom, key, value)
        # bpy.context.view_layer.objects.active = batom.obj

def get_scale(self):
    batom = get_active_batom()
    if batom is not None:
        return batom.scale[0]
    else:
        return 0

def set_scale(self, value):
    self["scale"] = value
    set_batom_attr("scale", value)




class BatomProperties(bpy.types.PropertyGroup):

    def Callback_modify_size(self, context):
        batom = bpy.context.scene.batoms.batom
        size = batom.size
        bpy.ops.batoms.batom_modify(key='size', size=size)

    def Callback_modify_bond(self, context):
        batom = bpy.context.scene.batoms.batom
        bond = batom.bond
        bpy.ops.batoms.batom_modify(key='bond', bond=bond)

    def Callback_modify_show(self, context):
        batom = bpy.context.scene.batoms.batom
        show = batom.show
        bpy.ops.batoms.batom_modify(key='show', show=show)

    scale: FloatProperty(
        name="scale", default=0.6,
        min=0, soft_max=2,
        description="scale",
        get=get_scale,
        set=set_scale)

    size: FloatProperty(
        name="size", default=1.5,
        min=0, soft_max=4,
        description="size", update=Callback_modify_size)

    bond: BoolProperty(name="Bond", default=True,
                       description="bond",
                       update=Callback_modify_bond)

    show: BoolProperty(name="Show", default=True,
                       description="show",
                       update=Callback_modify_show)
