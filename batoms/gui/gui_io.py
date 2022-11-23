"""
"""
import bpy
from bpy.types import Panel, Operator
from bpy.props import StringProperty

class BatomsPropertiesIO(bpy.types.PropertyGroup):
    pubchem_cid: StringProperty(
        name="pubchem_cid", default='31423',
        description="pubchem_cid")

    mp_key: StringProperty(
        name="mp_Key", default='',
        description="mp_key")

    mp_id: StringProperty(
        name="mp_ID", default='mp-2815',
        description="mp_id")

    rscb_name: StringProperty(
        name="rscb_name", default='1ema',
        description="rscb_name")

class Batoms_IO_search_material_project(Operator):
    bl_idname = "batoms.search_materials_project"
    bl_label = "Search"
    bl_description = ("Search structure by id")

    def execute(self, context):
        from batoms.database.pymatgen import pymatgen_search
        io = context.scene.batoms.io
        batoms = pymatgen_search(io.mp_key, io.mp_id)
        return {'FINISHED'}


class VIEW3D_PT_Batoms_io_materials_project(Panel):
    bl_label = "Materials Project"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "VIEW3D_PT_Batoms_io_materials_project"
    bl_parent_id = 'VIEW3D_PT_Batoms_io'

    def draw(self, context):
        layout = self.layout
        io = context.scene.batoms.io

        layout.label(text="Search structure")
        layout.prop(io, "mp_key")
        layout.prop(io, "mp_id")
        layout.operator("batoms.search_materials_project")


class Batoms_IO_search_pubchem(Operator):
    bl_idname = "batoms.search_pubchem"
    bl_label = "Search"
    bl_description = ("Search structure by id")

    def execute(self, context):
        from batoms.database.pubchem import pubchem_search
        io = context.scene.batoms.io
        batoms = pubchem_search(io.pubchem_cid)
        return {'FINISHED'}


class VIEW3D_PT_Batoms_io_pubchem(Panel):
    bl_label = "Pubchem"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "VIEW3D_PT_Batoms_io_pubchem"
    bl_parent_id = 'VIEW3D_PT_Batoms_io'

    def draw(self, context):
        layout = self.layout
        io = context.scene.batoms.io

        layout.label(text="Search structure")
        layout.prop(io, "pubchem_cid")
        layout.operator("batoms.search_pubchem")


class Batoms_IO_search_RSCB(Operator):
    bl_idname = "batoms.search_rscb"
    bl_label = "Import"
    bl_description = ("Import pdb by name")

    def execute(self, context):
        from batoms.database.rscb import rscb_import
        io = context.scene.batoms.io
        batoms = rscb_import(io.rscb_name)
        batoms.ribbon.draw()
        return {'FINISHED'}


class VIEW3D_PT_Batoms_io_RSCB(Panel):
    bl_label = "RSCB"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "VIEW3D_PT_Batoms_io_RSCB"
    bl_parent_id = 'VIEW3D_PT_Batoms_io'

    def draw(self, context):
        layout = self.layout
        io = context.scene.batoms.io

        layout.label(text="Import PDB")
        layout.prop(io, "rscb_name")
        layout.operator("batoms.search_rscb")

class VIEW3D_PT_Batoms_io(Panel):
    bl_label = "Import-Export"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "VIEW3D_PT_Batoms_io"

    def draw(self, context):
        layout = self.layout
        layout.label(text="Import")
        layout.operator("batoms.import")
        layout.label(text="Export")
        layout.operator("batoms.export")
