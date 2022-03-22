"""
this adds your menu to shift-a add object menu
please look at other menu classes in \scripts\startup\bl_ui\
i.e. space_view3d.py
"""
import bpy
from bpy.types import Menu


def menu_func(self, context):
    self.layout.menu("VIEW3D_MT_batoms_add", icon="MESH_UVSPHERE")
    self.layout.menu("VIEW3D_MT_surface_add", icon="MOD_LATTICE")
    self.layout.menu("VIEW3D_MT_nanotube_add", icon="META_CAPSULE")
    self.layout.menu("VIEW3D_MT_nanoparticle_add", icon="MESH_ICOSPHERE")


class VIEW3D_MT_batoms_add(Menu):
    bl_idname = "VIEW3D_MT_batoms_add"
    bl_label = "Batoms"
    bl_icon = "batoms_molecule"

    def draw(self, _context):
        layout = self.layout

        layout.operator_context = 'INVOKE_REGION_WIN'

        layout.operator("batoms.atoms_add", text="Atoms", icon='LAYER_ACTIVE')
        layout.operator("batoms.molecule_add", text="Molecule",
                        icon='FORCE_LENNARDJONES')
        layout.operator("batoms.bulk_add", text="Bulk", icon='MESH_CUBE')
        layout.operator("batoms.surface_add",
                        text="Surface", icon='VIEW_ORTHO')
        layout.operator("batoms.root_surface_add",
                        text="Root Surface", icon='MOD_LATTICE')
        layout.operator("batoms.smiles_add",
                        text="Smiles", icon='MOD_LATTICE')

        layout.separator()


class VIEW3D_MT_surface_add(Menu):
    bl_idname = "VIEW3D_MT_surface_add"
    bl_label = "Surface"
    # bl_icon = ""

    def draw(self, _context):
        layout = self.layout

        layout.operator_context = 'INVOKE_REGION_WIN'

        layout.operator("surface.fcc100_add", text="fcc100", icon='MESH_PLANE')
        layout.operator("surface.fcc110_add", text="fcc110", icon='MESH_PLANE')
        layout.operator("surface.fcc111_add", text="fcc111", icon='MESH_PLANE')
        layout.operator("surface.fcc211_add", text="fcc211", icon='MESH_PLANE')
        layout.operator("surface.fcc111_root_add",
                        text="fcc111_root", icon='MESH_PLANE')
        layout.separator()
        layout.operator("surface.bcc100_add", text="bcc100", icon='MESH_PLANE')
        layout.operator("surface.bcc110_add", text="bcc110", icon='MESH_PLANE')
        layout.operator("surface.bcc111_add", text="bcc111", icon='MESH_PLANE')
        layout.operator("surface.bcc111_root_add",
                        text="bcc111_root", icon='MESH_PLANE')
        layout.separator()
        layout.operator("surface.hcp0001_add",
                        text="hcp0001", icon='MESH_PLANE')
        layout.operator("surface.hcp10m10_add",
                        text="hcp10m10", icon='MESH_PLANE')
        layout.operator("surface.hcp0001_root_add",
                        text="hcp0001_root", icon='MESH_PLANE')
        layout.separator()
        layout.operator("surface.diamond100_add",
                        text="diamond100", icon='MESH_PLANE')
        layout.operator("surface.diamond111_add",
                        text="diamond111", icon='MESH_PLANE')


class VIEW3D_MT_nanotube_add(Menu):
    bl_idname = "VIEW3D_MT_nanotube_add"
    bl_label = "Nanotube"
    # bl_icon = ""

    def draw(self, _context):
        layout = self.layout

        layout.operator_context = 'INVOKE_REGION_WIN'

        layout.operator("nano.nanotube_add", text="Nanotube",
                        icon='MESH_CYLINDER')
        layout.operator("nano.nanoribbon_add",
                        text="nanoribbon", icon='MESH_PLANE')


class VIEW3D_MT_nanoparticle_add(Menu):
    bl_idname = "VIEW3D_MT_nanoparticle_add"
    bl_label = "Nanoparticle"
    # bl_icon = ""

    def draw(self, _context):
        layout = self.layout

        layout.operator_context = 'INVOKE_REGION_WIN'

        layout.operator("nano.decahedron_add",
                        text="Decahedron", icon='MESH_ICOSPHERE')
        layout.operator("nano.icosahedron_add",
                        text="Icosahedron", icon='MESH_ICOSPHERE')
        layout.operator("nano.octahedron_add",
                        text="Octahedron", icon='MESH_ICOSPHERE')
