import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (StringProperty,
                       )
from batoms import Batoms
from batoms.plugins.pymatgen import pymatgen_search

# The panel.
class Pymatgen_PT_prepare(Panel):
    bl_label       = "Pymatgen"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "PYMATGEN_PT_Tools"

    def draw(self, context):
        layout = self.layout
        pmgpanel = context.scene.pmgpanel

        layout.label(text="Search structure")
        layout.prop(pmgpanel, "key")
        layout.prop(pmgpanel, "id")
        layout.operator("batoms.pymatgen_search")


class PymatgenProperties(bpy.types.PropertyGroup):
    
    key: StringProperty(
        name = "Key", default='',
        description = "key")
    id: StringProperty(
        name = "ID", default='mp-2815',
        description = "id")

class Search(Operator):
    bl_idname = "batoms.pymatgen_search"
    bl_label = "Search"
    bl_description = ("Search structure by id")
    def execute(self, context):
        pmgpanel = context.scene.pmgpanel
        batoms = pymatgen_search(pmgpanel.key, pmgpanel.id)
        return {'FINISHED'}

