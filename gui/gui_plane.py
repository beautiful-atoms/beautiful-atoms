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
class Plane_PT_prepare(Panel):
    bl_label       = "Plane Tools"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    # bl_options     = {}
    bl_category = "Plane"
    bl_idname = "PLANE_PT_Tools"

    def draw(self, context):
        layout = self.layout
        plpanel = context.scene.plpanel

        box = layout.box()
        row = box.row()
        row.prop(plpanel, "indices")
        row = box.row()
        row.prop(plpanel, "distance")
        col = box.column(align=True)
        row = box.row()
        row.prop(plpanel, "color")
        box = layout.box()
        col = box.column(align=True)
        col.operator("batoms.add_plane")



class PlaneProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    @property
    def selected_plane(self):
        return get_selected_objects('bplane')
    def Callback_modify_indices(self, context):
        plpanel = bpy.context.scene.plpanel
        modify_plane_attr(self.selected_batoms, self.selected_plane, 'indices', plpanel.indices)
    def Callback_modify_distance(self, context):
        plpanel = bpy.context.scene.plpanel
        modify_plane_attr(self.selected_batoms, self.selected_plane, 'distance', plpanel.distance)
    def Callback_modify_color(self, context):
        plpanel = bpy.context.scene.plpanel
        color = plpanel.color
        modify_plane_attr(self.selected_batoms, self.selected_plane, 'color', color)

    
    indices: IntVectorProperty(
        name="Miller indices", default=(0, 0, 1),
        description = "Miller indices for the plane", update = Callback_modify_indices)
    distance: FloatProperty(
        name="distance", default=0.5,
        description = "distance from origin", update = Callback_modify_distance)
    color: FloatVectorProperty(
        name="color", 
        subtype='COLOR',
        default=(0.1, 0.8, 0.4, 0.8),
        size =4,
        description="color picker",
        update = Callback_modify_color)
    

def modify_plane_attr(batoms_name_list, plpanel_name_list, key, value):
    selected_plane_new = []
    for batoms_name in batoms_name_list:
        batoms = Batoms(label = batoms_name)
        for plane_name in plpanel_name_list:
            plane = bpy.data.objects[plane_name]
            if plane.bplane.label == batoms_name:
                setattr(batoms.planesetting[plane.bplane.indices], key, value)
            selected_plane_new.append(plane_name)
        batoms.draw_lattice_plane()
    for name in selected_plane_new:
        obj = bpy.data.objects.get(name)
        if obj is not None:
            obj.select_set(True)        
    

def add_plane(indices, color, distance):
    """
    """
    selected_batoms = get_selected_batoms()
    for batoms_name in selected_batoms:
        batoms = Batoms(label = batoms_name)
        batoms.planesetting[indices] = {'indices': indices, 
                                'distance': distance,
                                'color': color}
        batoms.draw_lattice_plane()

class AddButton(Operator):
    bl_idname = "batoms.add_plane"
    bl_label = "Add"
    bl_description = "Add distance, angle and dihedra angle"

    def execute(self, context):
        plpanel = context.scene.plpanel
        add_plane(plpanel.indices, plpanel.color, plpanel.distance)
        return {'FINISHED'}