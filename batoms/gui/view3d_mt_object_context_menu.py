"""
view3d_mt_object_context_menu
"""

import bpy
from bpy.types import Menu


def menu_func(self, context):
    self.layout.menu("VIEW3D_MT_object_context_batoms", icon="MESH_UVSPHERE")


class VIEW3D_MT_object_context_batoms(Menu):
    bl_idname = "VIEW3D_MT_object_context_batoms"
    bl_label = "Batoms"
    bl_icon = "batoms_molecule"

    def draw(self, _context):
        layout = self.layout

        layout.operator_context = 'INVOKE_REGION_WIN'

        layout.operator("batoms.join", text="Join", icon='LAYER_ACTIVE')
        layout.operator("batoms.delete_selected", text="Delete",
                        icon='FORCE_LENNARDJONES')

        layout.separator()