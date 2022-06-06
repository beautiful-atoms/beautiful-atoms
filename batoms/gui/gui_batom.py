import bpy
from bpy.types import (Panel,
                       )
from bpy.props import (FloatProperty,
                       FloatVectorProperty,
                       BoolProperty,
                       EnumProperty,
                       )
from batoms.gui.gui_batoms import get_active_batoms
from batoms.utils.butils import get_selected_vertices_bmesh

class Batom_PT_prepare(Panel):
    bl_label = "Batom"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Species"
    bl_idname = "BATOMS_PT_Batom"

    def draw(self, context):
        layout = self.layout
        btpanel = context.scene.btpanel
        layout.prop(btpanel, "scale")
        layout.prop(btpanel, "show")
        layout.prop(btpanel, "bond")



def get_scale(self):
    batoms = get_active_batoms()
    if batoms is not None:
        v = get_selected_vertices_bmesh(batoms.obj)
        if len(v) > 0:
            return batoms.scale[v[0]]
    else:
        return 0

def set_scale(self, value):
    self["scale"] = value
    bpy.ops.batoms.batom_modify(key='scale', scale=value)

class BatomProperties(bpy.types.PropertyGroup):

    def Callback_modify_size(self, context):
        btpanel = bpy.context.scene.btpanel
        size = btpanel.size
        bpy.ops.batoms.batom_modify(key='size', size=size)

    def Callback_modify_bond(self, context):
        btpanel = bpy.context.scene.btpanel
        bond = btpanel.bond
        bpy.ops.batoms.batom_modify(key='bond', bond=bond)

    def Callback_modify_show(self, context):
        btpanel = bpy.context.scene.btpanel
        show = btpanel.show
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
