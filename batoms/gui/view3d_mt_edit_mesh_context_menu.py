"""
VIEW3D_MT_edit_mesh_context_menu
"""

import bpy
from bpy.types import Menu


def menu_func(self, context):
    self.layout.menu("VIEW3D_MT_edit_mesh_context_batoms", icon="MESH_UVSPHERE")



class VIEW3D_MT_edit_mesh_context_batoms_model_style(Menu):
    bl_idname = "VIEW3D_MT_edit_mesh_context_batoms_model_style"
    bl_label = "model style"
    # bl_icon = ""

    def draw(self, context):
        is_vert_mode, is_edge_mode, is_face_mode = context.tool_settings.mesh_select_mode

        layout = self.layout

        layout.operator_context = 'INVOKE_REGION_WIN'

        op0 = layout.operator("batoms.apply_model_style_selected", text="Space-filling")
        op0.model_style = '0'
        op1 = layout.operator("batoms.apply_model_style_selected", text="Ball-and-stick")
        op1.model_style = '1'
        op2 = layout.operator("batoms.apply_model_style_selected", text="Polyhedral")
        op2.model_style = '2'
        op3 = layout.operator("batoms.apply_model_style_selected", text="Stick")
        op3.model_style = '3'



class VIEW3D_MT_edit_mesh_context_batoms(Menu):
    bl_idname = "VIEW3D_MT_edit_mesh_context_batoms"
    bl_label = "Batoms"
    bl_icon = "batoms_molecule"

    def draw(self, _context):
        layout = self.layout

        layout.operator_context = 'INVOKE_REGION_WIN'
        layout.menu("VIEW3D_MT_edit_mesh_context_batoms_model_style",
                    text="Model Style", icon='LAYER_ACTIVE')
        # layout.menu("VIEW3D_MT_edit_mesh_context_batoms_color_style",
                    # text="Color Style", icon='LAYER_ACTIVE')
        # layout.menu("VIEW3D_MT_edit_mesh_context_batoms_radius_style",
                    # text="Radius Style", icon='LAYER_ACTIVE')
        layout.operator("batoms.replace", text="Replace", icon='LAYER_ACTIVE')
        # layout.operator("batoms.select", text="Select", icon='LAYER_ACTIVE')
        layout.operator("batoms.delete_selected_atoms", text="Delete",
                        icon='FORCE_LENNARDJONES')

        layout.separator()
