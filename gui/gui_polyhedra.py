import bpy
from bpy.types import (Panel,
                       )
from bpy.props import (BoolProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       )

from batoms.butils import get_selected_objects, get_selected_batoms
from batoms import Batoms
from batoms.gui.gui_batoms import modify_batoms_attr

class Polyhedra_PT_prepare(Panel):
    bl_label       = "Polyhedra"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "BATOMS_PT_Polyhedra"

  
    def draw(self, context):
        layout = self.layout
        popanel = context.scene.popanel
    
        layout.label(text="Polyhedra style")
        layout.prop(popanel, "polyhedra_style", expand  = True)

        layout.prop(popanel, "show_edge")
        layout.prop(popanel, "width")

        layout.prop(popanel, "polyhedracolor")

class PolyhedraProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    @property
    def selected_polyhedra(self):
        return get_selected_objects('bpolyhedra')
    def Callback_polyhedra_style(self, context):
        popanel = bpy.context.scene.popanel
        polyhedra_style = list(popanel.polyhedra_style)[0]
        modify_batoms_attr(self.selected_batoms, 'polyhedra_style', polyhedra_style)
    def Callback_modify_width(self, context):
        popanel = bpy.context.scene.popanel
        width = popanel.width
        modify_polyhedra_attr(self.selected_batoms, self.selected_polyhedra, 'width', width)
    def Callback_modify_show_edge(self, context):
        popanel = bpy.context.scene.popanel
        show_edge = popanel.show_edge
        modify_polyhedra_attr(self.selected_batoms, self.selected_polyhedra, 'show_edge', show_edge)
    def Callback_modify_polyhedracolor(self, context):
        popanel = bpy.context.scene.popanel
        polyhedracolor = popanel.polyhedracolor
        modify_polyhedra_attr(self.selected_batoms, self.selected_polyhedra, 'color', polyhedracolor)

    polyhedra_style: EnumProperty(
        name="polyhedra_style",
        description="Polhhedra models",
        items=(('0',"0", "atoms, bonds and polyhedra"),
               ('1',"1", "atoms, polyhedra"),
               ('2',"2","central atoms, polyhedra"),
               ('3',"3", "polyhedra")),
        default={'0'}, 
        update=Callback_polyhedra_style,
        options={'ENUM_FLAG'},
        )
    width: FloatProperty(
        name="edgewidth", default=0.01,
        description = "width", update = Callback_modify_width)
    show_edge: BoolProperty(
        name = "show_edge", default=True,
        description = "show_edge", update = Callback_modify_show_edge)
    polyhedracolor: FloatVectorProperty(
        name="polyhedracolor", 
        subtype='COLOR',
        default=(0.1, 0.8, 0.4 ,1.0),
        size =4,
        description="color picker",
        update = Callback_modify_polyhedracolor)

def modify_polyhedra_attr(selected_batoms, selected_polyhedra, key, value):
    selected_polyhedra_new = []
    for batoms_name in selected_batoms:
        batoms = Batoms(label = batoms_name)
        for polyhedra_name in selected_polyhedra:
            polyhedra = bpy.data.objects[polyhedra_name]
            if polyhedra.bpolyhedra.label == batoms_name:
                setattr(batoms.polyhedrasetting[polyhedra.bpolyhedra.species], key, value)
                selected_polyhedra_new.append(polyhedra_name)
        batoms.draw_polyhedras()
    for name in selected_polyhedra_new:
        obj = bpy.data.objects.get(name)
        obj.select_set(True)        