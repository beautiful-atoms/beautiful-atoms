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
    self.layout.menu("VIEW3D_MT_surface_add")
    self.layout.menu("VIEW3D_MT_nanotube_add")
    self.layout.menu("VIEW3D_MT_nanoparticle_add")


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


class VIEW3D_MT_surface_add(Menu):
    bl_idname = "VIEW3D_MT_surface_add"
    bl_label = "Surface"
    # bl_icon = ""

    def draw(self, _context):
        layout = self.layout

        layout.operator_context = 'INVOKE_REGION_WIN'

        layout.operator("surface.fcc100", text="FCC100", icon='MESH_PLANE')
        layout.operator("surface.fcc110", text="FCC110", icon='MESH_PLANE')
        layout.operator("surface.fcc111", text="FCC111", icon='MESH_PLANE')



class VIEW3D_MT_nanotube_add(Menu):
    bl_idname = "VIEW3D_MT_nanotube_add"
    bl_label = "Nanotube"
    # bl_icon = ""

    def draw(self, _context):
        layout = self.layout

        layout.operator_context = 'INVOKE_REGION_WIN'

        layout.operator("nanotube.nanotube", text="Nanotube", icon='MESH_CYLINDER')
        # layout.operator("nanotube.fcc110", text="FCC110", icon='MESH_PLANE')
        # layout.operator("nanotube.fcc111", text="FCC111", icon='MESH_PLANE')



class VIEW3D_MT_nanoparticle_add(Menu):
    bl_idname = "VIEW3D_MT_nanoparticle_add"
    bl_label = "Nanoparticle"
    # bl_icon = ""

    def draw(self, _context):
        layout = self.layout

        layout.operator_context = 'INVOKE_REGION_WIN'

        layout.operator("nanoparticle.decahedron", text="Decahedron", icon='MESH_CUBE')
        layout.operator("nanoparticle.icosahedron", text="Icosahedron", icon='MESH_CUBE')
        layout.operator("nanoparticle.octahedron", text="Octahedron", icon='MESH_CUBE')


