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

from batoms.butils import get_selected_objects, get_selected_batoms
from batoms import Batoms

class Polyhedra_PT_prepare(Panel):
    bl_label       = "Polyhedra"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    # bl_options     = {'DEFAULT_CLOSED'}
    bl_category = "Bond"
    bl_idname = "BATOMS_PT_Polyhedra"

  
    def draw(self, context):
        layout = self.layout
        plpanel = context.scene.plpanel

        
        box = layout.box()
        row = box.row()
        row.prop(plpanel, "edgewidth")

        col = box.column(align=True)
        row = box.row()
        row.prop(plpanel, "polyhedracolor")

class PolyhedraProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    @property
    def selected_polyhedra(self):
        return get_selected_objects('bpolyhedra')
    def Callback_modify_edgewidth(self, context):
        plpanel = bpy.context.scene.plpanel
        edgewidth = plpanel.edgewidth
        modify_polyhedra_attr(self.selected_batoms, self.selected_polyhedra, 'edgewidth', edgewidth)
    def Callback_modify_polyhedracolor(self, context):
        plpanel = bpy.context.scene.plpanel
        polyhedracolor = plpanel.polyhedracolor
        modify_polyhedra_attr(self.selected_batoms, self.selected_polyhedra, 'color', polyhedracolor)

    
    edgewidth: FloatProperty(
        name="Edgewidth", default=0.2,
        description = "edgewidth", update = Callback_modify_edgewidth)
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