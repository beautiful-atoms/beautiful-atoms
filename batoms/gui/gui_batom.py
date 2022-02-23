import bpy
from bpy.types import (Panel,
                       )
from bpy.props import (FloatProperty,
                       FloatVectorProperty,
                       StringProperty,
                       EnumProperty,
                       )

from batoms.utils.butils import get_selected_batoms, get_selected_vertices
from batoms import Batoms

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

        layout.label(text="Shape")
        col = layout.column()
        col.prop(btpanel, "batom_shape", expand  = True)
        layout.prop(btpanel, "elements")
        layout.prop(btpanel, "scale")
        layout.prop(btpanel, "size")

        layout.prop(btpanel, "batomcolor")

class BatomProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms('batom')
    @property
    def selected_vertices(self):
        return get_selected_vertices()
    def Callback_batom_shape(self, context):
        btpanel = bpy.context.scene.btpanel
        batom_shape = int(list(btpanel.batom_shape)[0])
        modify_batom_attr(self.selected_vertices, 'shape', batom_shape)
    def Callback_modify_scale(self, context):
        btpanel = bpy.context.scene.btpanel
        scale = btpanel.scale
        modify_batom_attr(self.selected_vertices, 'scale', scale)
    def Callback_modify_size(self, context):
        btpanel = bpy.context.scene.btpanel
        size = btpanel.size
        modify_batom_attr(self.selected_vertices, 'size', size)
    def Callback_modify_elements(self, context):
        import json
        btpanel = bpy.context.scene.btpanel
        elements = btpanel.elements
        elements = elements.replace("'", '"')
        elements = json.loads(elements)
        modify_batom_attr(self.selected_vertices, 'elements', elements)
    def Callback_modify_batomcolor(self, context):
        btpanel = bpy.context.scene.btpanel
        batomcolor = btpanel.batomcolor
        modify_batom_attr(self.selected_vertices, 'color', batomcolor)

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
        name="scale", default=0.6,
        description = "scale", update = Callback_modify_scale)
    size: FloatProperty(
        name="size", default=1.5,
        description = "size", update = Callback_modify_size)
    elements: StringProperty(
        name = "Elements", default='{"Al": 0.7, "Si": 0.3}',
        description = "str or dict", update = Callback_modify_elements)
    batomcolor: FloatVectorProperty(
        name="color", 
        subtype='COLOR',
        default=(0.1, 0.8, 0.4 ,1.0),
        size =4,
        description="color picker",
        update = Callback_modify_batomcolor)

def modify_batom_attr(selected_vertices, key, value):
    for label, indices in selected_vertices:
        batoms = Batoms(label = label)
        if key in ['scale', 'size']:
            array0 = batoms.attributes[key]
            array0[indices] = value
            batoms.set_attributes({key: array0})