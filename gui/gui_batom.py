import bpy
from bpy.types import (Panel,
                       Operator,
                       AddonPreferences,
                       PropertyGroup,
                       )
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       PointerProperty,
                       )

from batoms.butils import get_selected_objects
from batoms import Batom

class Batom_PT_prepare(Panel):
    bl_label       = "Batom"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "BATOMS_PT_Batom"

  
    def draw(self, context):
        layout = self.layout
        btpanel = context.scene.btpanel

        box = layout.box()
        col = box.column()
        col.label(text="Batom shape")
        col = box.column()
        col.prop(btpanel, "batom_shape", expand  = True)
        box = layout.box()
        col = box.column(align=True)
        row = box.row()
        row.prop(btpanel, "scale")

        col = box.column(align=True)
        row = box.row()
        row.prop(btpanel, "batomcolor")

class BatomProperties(bpy.types.PropertyGroup):
    @property
    def selected_batom(self):
        return get_selected_objects('batom')
    def Callback_batom_shape(self, context):
        btpanel = bpy.context.scene.btpanel
        batom_shape = int(list(btpanel.batom_shape)[0])
        modify_batom_attr(self.selected_batom, 'shape', batom_shape)
    def Callback_modify_scale(self, context):
        btpanel = bpy.context.scene.btpanel
        scale = btpanel.scale
        modify_batom_attr(self.selected_batom, 'scale', scale)
    def Callback_modify_batomcolor(self, context):
        btpanel = bpy.context.scene.btpanel
        batomcolor = btpanel.batomcolor
        modify_batom_attr(self.selected_batom, 'color', batomcolor)

    batom_shape: EnumProperty(
        name="shape",
        description="batom shape",
        items=(('0',"UV_Sphere", ""),
               ('1',"ICO_Sphere", ""),
               ('2',"Cube", "")),
        default={'0'},
        update=Callback_batom_shape,
        options={'ENUM_FLAG'},
        )
    scale: FloatProperty(
        name="scale", default=0.2,
        description = "scale", update = Callback_modify_scale)
    batomcolor: FloatVectorProperty(
        name="color", 
        subtype='COLOR',
        default=(0.1, 0.8, 0.4 ,1.0),
        size =4,
        description="color picker",
        update = Callback_modify_batomcolor)

def modify_batom_attr(batom_name_list, key, value):
    batom_list = []
    for name in batom_name_list:
        batom = Batom(label = name)
        setattr(batom, key, value)
        batom_list.append(batom)
    for batom in batom_list:
        batom.select = True