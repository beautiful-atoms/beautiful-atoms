import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (BoolProperty,
                       IntVectorProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       )
from batoms.utils.butils import get_selected_batoms, get_selected_objects
from batoms.batoms import Batoms

# The panel.
class Plane_PT_prepare(Panel):
    bl_label       = "Plane"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Surface"
    bl_idname = "PLANE_PT_Tools"

    def draw(self, context):
        plpanel = context.scene.plpanel
        
        layout = self.layout
        row = layout.row()
        row.prop(plpanel, "indices")
        layout.prop(plpanel, "distance")
        layout.prop(plpanel, "scale")
        layout.prop(plpanel, "symmetry", expand  = True)
        row = layout.row()
        row.prop(plpanel, "crystal")
        row.prop(plpanel, "center")
        layout.prop(plpanel, "slicing", expand  = True)
        layout.prop(plpanel, "boundary", expand  = True)
        layout.prop(plpanel, "show_edge")
        layout.prop(plpanel, "color")
        layout.operator("batoms.add_plane")



class PlaneProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    @property
    def selected_plane(self):
        return get_selected_objects('bplane')
    def Callback_modify_indices(self, context):
        plpanel = bpy.context.scene.plpanel
        modify_plane_attr(self.selected_batoms, 
                    self.selected_plane, 'indices', plpanel.indices)
    def Callback_modify_distance(self, context):
        plpanel = bpy.context.scene.plpanel
        modify_plane_attr(self.selected_batoms, 
                    self.selected_plane, 'distance', plpanel.distance)
    def Callback_modify_scale(self, context):
        plpanel = bpy.context.scene.plpanel
        modify_plane_attr(self.selected_batoms, 
                    self.selected_plane, 'scale', plpanel.scale)
    def Callback_modify_crystal(self, context):
        plpanel = bpy.context.scene.plpanel
        modify_plane_attr(self.selected_batoms, 
                    self.selected_plane, 'crystal', plpanel.crystal, center = plpanel.center)
    def Callback_modify_center(self, context):
        plpanel = bpy.context.scene.plpanel
        modify_plane_attr(self.selected_batoms, 
                    self.selected_plane, 'center', plpanel.center)
    def Callback_modify_symmetry(self, context):
        plpanel = bpy.context.scene.plpanel
        modify_plane_attr(self.selected_batoms, 
                    self.selected_plane, 'symmetry', plpanel.symmetry)
    def Callback_modify_slicing(self, context):
        plpanel = bpy.context.scene.plpanel
        modify_plane_attr(self.selected_batoms, 
                    self.selected_plane, 'slicing', plpanel.slicing)
    def Callback_modify_boundary(self, context):
        plpanel = bpy.context.scene.plpanel
        modify_plane_attr(self.selected_batoms, 
                    self.selected_plane, 'boundary', plpanel.boundary)
    def Callback_modify_show_edge(self, context):
        plpanel = bpy.context.scene.plpanel
        show_edge = plpanel.show_edge
        modify_plane_attr(self.selected_batoms, 
                    self.selected_plane, 'show_edge', show_edge)
    def Callback_modify_color(self, context):
        plpanel = bpy.context.scene.plpanel
        color = plpanel.color
        modify_plane_attr(self.selected_batoms, 
                    self.selected_plane, 'color', color)

    
    indices: IntVectorProperty(
        name="Miller indices", size = 3, default=(0, 0, 1),
        description = "Miller indices for the plane", update = Callback_modify_indices)
    distance: FloatProperty(
        name="Distance", default=3,
        description = "distance from origin", update = Callback_modify_distance)
    scale: FloatProperty(
        name="Scale", default=1,
        description = "scale of the plane", update = Callback_modify_scale)
    crystal: BoolProperty(name="crystal",
                default=False, 
                description = "plane to form crystal shape",
                update = Callback_modify_crystal)
    center: BoolProperty(name="center",
                default=False, 
                description = "Apply center to model",
                update = Callback_modify_center)
    symmetry: BoolProperty(name="symmetry",
                default=False, 
                description = "Apply symmetry to indices",
                update = Callback_modify_symmetry)
    slicing: BoolProperty(name="slicing",
                default=False, 
                description = "Apply slicing to volumetric data",
                update = Callback_modify_slicing)
    boundary: BoolProperty(name="boundary",
                default=False, 
                description = "Apply boundary to model",
                update = Callback_modify_boundary)
    show_edge: BoolProperty(
        name = "show_edge", default=False,
        description = "show_edge", update = Callback_modify_show_edge)
    color: FloatVectorProperty(
        name="color", 
        subtype='COLOR',
        default=(0.1, 0.8, 0.4, 0.8),
        size =4,
        description="color picker",
        update = Callback_modify_color)
    

def modify_plane_attr(batoms_name_list, plpanel_name_list, key, value, center = False):
    selected_plane_new = []
    for batoms_name in batoms_name_list:
        batoms = Batoms(label = batoms_name)
        for plane_name in plpanel_name_list:
            plane = bpy.data.objects[plane_name]
            if plane.bplane.label == batoms_name:
                batoms.planesetting[plane.bplane.indices] = {key: value}
            selected_plane_new.append(plane_name)
        batoms.planesetting.draw_lattice_plane()
        if center:
            batoms.planesetting.draw_crystal_shape(origin = batoms.cell.center)
        else:
            batoms.planesetting.draw_crystal_shape()
    for name in selected_plane_new:
        obj = bpy.data.objects.get(name)
        if obj is not None:
            obj.select_set(True)      
    

def add_plane(indices, color, distance, scale, crystal, 
                symmetry, slicing, boundary, show_edge,
                center = False):
    """
    """
    selected_batoms = get_selected_batoms()
    for batoms_name in selected_batoms:
        batoms = Batoms(label = batoms_name)
        batoms.planesetting[indices] = {'indices': indices, 
                                'distance': distance,
                                'scale': scale,
                                'crystal': crystal,
                                'symmetry': symmetry,
                                'slicing': slicing,
                                'boundary': boundary,
                                'show_edge': show_edge,
                                'color': color, 
                                }
        batoms.planesetting.draw_lattice_plane()
        if center:
            batoms.planesetting.draw_crystal_shape(origin = batoms.cell.center)
        else:
            batoms.planesetting.draw_crystal_shape()

class AddButton(Operator):
    bl_idname = "batoms.add_plane"
    bl_label = "Add"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Add plane"

    indices: IntVectorProperty(
        name="Miller indices", size = 3, default=(0, 0, 1),
        description = "Miller indices for the plane", )
    distance: FloatProperty(
        name="Distance", default=3,
        description = "distance from origin", )
    scale: FloatProperty(
        name="Scale", default=1,
        description = "scale of the plane", )
    crystal: BoolProperty(name="crystal",
                default=False, 
                description = "plane to form crystal shape",
                )
    center: BoolProperty(name="center",
                default=False, 
                description = "Apply center to model",
                )
    symmetry: BoolProperty(name="symmetry",
                default=False, 
                description = "Apply symmetry to indices",
                )
    slicing: BoolProperty(name="slicing",
                default=False, 
                description = "Apply slicing to volumetric data",
                )
    boundary: BoolProperty(name="boundary",
                default=False, 
                description = "Apply boundary to model",
                )
    show_edge: BoolProperty(
        name = "show_edge", default=False,
        description = "show_edge", )
    color: FloatVectorProperty(
        name="color", 
        subtype='COLOR',
        default=(0.1, 0.8, 0.4, 0.8),
        size =4,
        description="color picker",
        )

    def execute(self, context):
        add_plane(self.indices, 
                    self.color, 
                    self.distance, 
                    self.scale, 
                    self.crystal,
                    self.symmetry,
                    self.slicing,
                    self.boundary,
                    self.show_edge,
                    self.center,
                    )
        return {'FINISHED'}