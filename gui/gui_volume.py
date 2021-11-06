import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (StringProperty,
                       BoolProperty,
                       BoolVectorProperty,
                       IntProperty,
                       IntVectorProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       PointerProperty,
                       )
from batoms.butils import get_selected_batoms, get_selected_objects
from batoms.batoms import Batoms

# The panel.
class Volume_PT_prepare(Panel):
    bl_label       = "Volume"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "VOLUME_PT_Tools"

    def draw(self, context):
        layout = self.layout
        vopanel = context.scene.vopanel

        box = layout.box()
        row = box.row()
        row.prop(vopanel, "level")
        col = box.column(align=True)
        row = box.row()
        row.prop(vopanel, "color")

        box = layout.box()
        col = box.column(align=True)
        col.operator("batoms.add_level")
        col.prop(vopanel, "label")



class VolumeProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    @property
    def selected_isosurface(self):
        return get_selected_objects('bisosurface')
    def Callback_modify_level(self, context):
        vopanel = bpy.context.scene.vopanel
        modify_volume_attr(self.selected_batoms, self.selected_isosurface, 'level', vopanel.level)
    def Callback_modify_color(self, context):
        vopanel = bpy.context.scene.vopanel
        color = vopanel.color
        modify_volume_attr(self.selected_batoms, self.selected_isosurface, 'color', color)

    
    level: FloatProperty(
        name="Level", default=0.01,
        description = "level for isosurface", update = Callback_modify_level)
    color: FloatVectorProperty(
        name="color", 
        subtype='COLOR',
        default=(0.1, 0.8, 0.4 ,1.0),
        size =4,
        description="color picker",
        update = Callback_modify_color)
    label: StringProperty(
        name="label", default='2',
        description = "new level")

def modify_volume_attr(batoms_name_list, bisosurface_name_list, key, value):
    selected_isosurface_new = []
    for batoms_name in batoms_name_list:
        batoms = Batoms(label = batoms_name)
        for isosurface_name in bisosurface_name_list:
            iso = bpy.data.objects[isosurface_name]
            if iso.bisosurface.label == batoms_name:
                setattr(batoms.isosurfacesetting[iso.bisosurface.name], key, value)
            selected_isosurface_new.append(isosurface_name)
        batoms.draw_isosurface()
    for name in selected_isosurface_new:
        obj = bpy.data.objects.get(name)
        obj.select_set(True)        
    

def add_level(name, level, color):
    """
    """
    selected_batoms = get_selected_batoms()
    for batoms_name in selected_batoms:
        batoms = Batoms(label = batoms_name)
        batoms.isosurfacesetting[name] = {'level': level, 'color': color}
        batoms.draw_isosurface()

class AddButton(Operator):
    bl_idname = "batoms.add_level"
    bl_label = "Add"
    bl_description = "Add distance, angle and dihedra angle"

    def execute(self, context):
        vopanel = context.scene.vopanel
        add_level(vopanel.label, vopanel.level, vopanel.color)
        return {'FINISHED'}