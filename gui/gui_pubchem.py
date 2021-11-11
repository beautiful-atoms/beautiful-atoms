import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (StringProperty,
                       )
from batoms.plugins.pubchem import pubchem_search

# The panel.
class Pubchem_PT_prepare(Panel):
    bl_label       = "Pubchem"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "PUBCHEM_PT_Tools"

    def draw(self, context):
        layout = self.layout
        pubcpanel = context.scene.pubcpanel

        box = layout.box()
        col = box.column(align=True)
        col.label(text="Search structure")
        col.prop(pubcpanel, "cid")
        col.operator("batoms.pubchem_search")


class PubchemProperties(bpy.types.PropertyGroup):
    
    cid: StringProperty(
        name = "cid", default='31423',
        description = "cid")


class Search(Operator):
    bl_idname = "batoms.pubchem_search"
    bl_label = "Search"
    bl_description = ("Search structure by id")
    def execute(self, context):
        pubcpanel = context.scene.pubcpanel
        batoms = pubchem_search(pubcpanel.cid)
        return {'FINISHED'}

