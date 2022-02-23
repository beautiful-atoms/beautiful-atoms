import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (StringProperty,
                       )
from batoms.plugins.rscb import rscb_import

# The panel.
class RSCB_PT_prepare(Panel):
    bl_label       = "RSCB"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "RSCB_PT_Tools"

    def draw(self, context):
        layout = self.layout
        pubcpanel = context.scene.pubcpanel

        layout.label(text="Import PDB")
        layout.prop(pubcpanel, "name")
        layout.operator("batoms.rscb_import")


class RSCBProperties(bpy.types.PropertyGroup):
    
    name: StringProperty(
        name = "name", default='1ema',
        description = "name")

class RSCB_Import(Operator):
    bl_idname = "batoms.rscb_import"
    bl_label = "Import"
    bl_description = ("Import pdb by name")
    def execute(self, context):
        pubcpanel = context.scene.pubcpanel
        batoms = rscb_import(pubcpanel.name)
        batoms.ribbon.draw()
        return {'FINISHED'}

