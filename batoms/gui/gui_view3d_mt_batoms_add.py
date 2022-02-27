"""
this adds your menu to shift-a add object menu
please look at other menu classes in \scripts\startup\bl_ui\
i.e. space_view3d.py
"""
import bpy

from bpy.types import (
    Menu,
)


def menu_func(self, context):
    self.layout.menu("VIEW3D_MT_batoms_add")


class VIEW3D_MT_batoms_add(Menu):
    bl_idname = "VIEW3D_MT_batoms_add"
    bl_label = "Batoms"
    bl_icon = "batoms_add_molecule"

    def draw(self, _context):
        layout = self.layout

        layout.operator_context = 'INVOKE_REGION_WIN'

        layout.operator("batoms.add_atoms", text="Atoms", icon='MESH_PLANE')
        layout.operator("batoms.add_molecule", text="Molecule", icon='MESH_CUBE')
        layout.operator("batoms.add_bulk", text="Bulk", icon='MESH_CIRCLE')
        layout.operator("batoms.add_surface", text="Surface", icon='MESH_UVSPHERE')
        
        layout.separator()

        # layout.operator("batoms.add_surface", text="Grid", icon='MESH_GRID')

