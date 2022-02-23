import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (StringProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       )
from batoms.utils.butils import get_selected_batoms, get_selected_objects
from batoms.batoms import Batoms

# The panel.
class Volume_PT_prepare(Panel):
    bl_label       = "Surface"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "VOLUME_PT_Tools"

    def draw(self, context):
        layout = self.layout
        vopanel = context.scene.vopanel

        layout.label(text="Isosurface")
        layout.prop(vopanel, "level")
        layout.prop(vopanel, "color_iso")

        layout.prop(vopanel, "label")
        layout.operator("batoms.add_level")

        layout.label(text="SAS")
        layout.prop(vopanel, "color_sas")




class VolumeProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    @property
    def selected_isosurface(self):
        return get_selected_objects('bisosurface')
    def Callback_modify_level(self, context):
        vopanel = bpy.context.scene.vopanel
        modify_isosurface_attr(self.selected_batoms, self.selected_isosurface, 'level', vopanel.level)
    def Callback_modify_color_iso(self, context):
        vopanel = bpy.context.scene.vopanel
        color = vopanel.color_iso
        modify_isosurface_attr(self.selected_batoms, self.selected_isosurface, 'color', color)
    def Callback_modify_color_sas(self, context):
        vopanel = bpy.context.scene.vopanel
        color = vopanel.color_sas
        modify_ms_attr(self.selected_batoms, self.selected_isosurface, 'color', color)

    
    level: FloatProperty(
        name="Level", default=0.01,
        description = "level for isosurface", update = Callback_modify_level)
    color_iso: FloatVectorProperty(
        name="color", 
        subtype='COLOR',
        default=(0.1, 0.8, 0.4 ,1.0),
        size =4,
        description="color picker",
        update = Callback_modify_color_iso)
    color_sas: FloatVectorProperty(
        name="color", 
        subtype='COLOR',
        default=(0.1, 0.8, 0.4 ,1.0),
        size =4,
        description="color picker",
        update = Callback_modify_color_sas)
    label: StringProperty(
        name="label", default='2',
        description = "new level")

def modify_isosurface_attr(batoms_name_list, bisosurface_name_list, key, value):
    selected_isosurface_new = []
    for batoms_name in batoms_name_list:
        batoms = Batoms(label = batoms_name)
        for isosurface_name in bisosurface_name_list:
            iso = bpy.data.objects[isosurface_name]
            if iso.batoms.bisosurface.label == batoms_name:
                setattr(batoms.isosurfacesetting[iso.batoms.bisosurface.name], key, value)
            selected_isosurface_new.append(isosurface_name)
        batoms.isosurfacesetting.draw_isosurface()
    for name in selected_isosurface_new:
        obj = bpy.data.objects.get(name)
        obj.select_set(True)

def modify_ms_attr(batoms_name_list, bms_name_list, key, value):
    selected_ms_new = []
    for batoms_name in batoms_name_list:
        batoms = Batoms(label = batoms_name)
        for ms_name in bms_name_list:
            obj = bpy.data.objects[ms_name]
            if obj.bms.label == batoms_name:
                setattr(batoms.mssetting[obj.bms.name], key, value)
            selected_ms_new.append(ms_name)
        batoms.mssetting.draw_SAS()

def add_level(name, level, color):
    """
    """
    selected_batoms = get_selected_batoms()
    for batoms_name in selected_batoms:
        batoms = Batoms(label = batoms_name)
        batoms.isosurfacesetting[name] = {'level': level, 'color': color}
        batoms.isosurfacesetting.draw_isosurface()

class AddButton(Operator):
    bl_idname = "batoms.add_level"
    bl_label = "Add"
    bl_description = "Add distance, angle and dihedra angle"

    def execute(self, context):
        vopanel = context.scene.vopanel
        add_level(vopanel.label, vopanel.level, vopanel.color)
        return {'FINISHED'}