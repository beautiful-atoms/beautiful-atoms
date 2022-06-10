import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (StringProperty,
                       )
from batoms.plugins.rscb import rscb_import


class RSCB_PT_prepare(Panel):
    bl_label = "RSCB"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "RSCB_PT_Tools"

    def draw(self, context):
        layout = self.layout
        rscb = context.scene.batoms.rscb

        layout.label(text="Import PDB")
        layout.prop(rscb, "name")
        layout.operator("batoms.rscb_import")


class RSCBProperties(bpy.types.PropertyGroup):

    name: StringProperty(
        name="name", default='1ema',
        description="name")


class RSCB_Import(Operator):
    bl_idname = "batoms.rscb_import"
    bl_label = "Import"
    bl_description = ("Import pdb by name")

    def execute(self, context):
        rscb = context.scene.batoms.rscb
        batoms = rscb_import(rscb.name)
        batoms.ribbon.draw()
        return {'FINISHED'}
