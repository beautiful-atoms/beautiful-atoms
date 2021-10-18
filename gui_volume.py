from ase.build import molecule, bulk
from ase import Atom, Atoms
import bpy
import bmesh
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
from batoms.butils import read_batoms_select
from batoms.batoms import Batoms

# The panel.
class Volume_PT_prepare(Panel):
    bl_label       = "Volume Tools"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    # bl_options     = {}
    bl_category = "Batoms"
    bl_idname = "VOLUME_PT_Tools"

    def draw(self, context):
        layout = self.layout
        blpanel = context.scene.blpanel

        box = layout.box()
        row = box.row()
        row.prop(blpanel, "level")



class VolumeProperties(bpy.types.PropertyGroup):
    @property
    def batoms_list(self):
        return self.get_batoms_list()
    def get_batoms_list(self):
        batoms_list = read_batoms_select()
        return batoms_list
    
    def Callback_modify_level(self, context):
        blpanel = bpy.context.scene.blpanel
        print('Callback_modify_level')
        modify_level(self.batoms_list, blpanel.level)
    
    level: FloatProperty(
        name="Level", default=0.01,
        description = "level for isosurface", update = Callback_modify_level)
    

def modify_level(batoms_name_list, level):
    for name in batoms_name_list:
        batoms = Batoms(label = name)
        batoms.isosurfacesetting.level = level
