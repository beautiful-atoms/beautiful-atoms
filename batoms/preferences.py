"""

0) update batoms
1) set custom file folder.
2) set custom windows

"""

import bpy
from bpy.types import AddonPreferences
from bpy.props import (
    BoolProperty,
)
from . import install_dependencies
import pathlib


cwd = pathlib.Path().resolve()

DEFAULT_GITHUB_ACCOUNT = "superstar54"
DEFAULT_REPO_NAME = "beautiful-atoms"
DEFAULT_PLUGIN_NAME = "batoms"

bpy.context.preferences.view.use_translate_new_dataname = False
bpy.context.preferences.inputs.use_rotate_around_active = True
bpy.context.preferences.inputs.use_zoom_to_mouse = True
# For laptop
bpy.context.preferences.inputs.use_emulate_numpad = True
# For laptop without mouse
bpy.context.preferences.inputs.use_mouse_emulate_3_button = True


# use default batoms interface
# set python console

# bpy.context.window.workspace = bpy.data.workspaces['UV Editing']





class BatomsAddonPreferences(AddonPreferences):
    bl_idname = __package__

    def Callback_install_ase(self, context):
        if self.has_ase:
            install_dependencies.install(["ase"])

    def Callback_install_scikit_image(self, context):
        if self.has_scikit_image:
            install_dependencies.install(["scikit-image"])

    def Callback_install_spglib(self, context):
        if self.has_spglib:
            install_dependencies.install(["spglib"])


    has_ase: BoolProperty(
        name="ASE installed",
        description="ASE package installed",
        default=False,
    )

    has_scikit_image: BoolProperty(
        name="scikit-image installed",
        description="scikit-image package installed",
        default=False,
    )

    has_spglib: BoolProperty(
        name="spglib installed",
        description="spglib package installed",
        default=False,
    )

    has_pymatgen: BoolProperty(
        name="pymatgen installed",
        description="pymatgen package installed",
        default=False,
    )

    
    def draw(self, context):
        layout = self.layout

        layout.label(text="Welcome to Batoms!")
        # Check Blender version
        if bpy.app.version_string < '3.0.0':
            box = layout.box().column()
            box.label(text="Warning: Batoms need Blender version > 3.0.0.")

        box = layout.box().column()
        row = box.row(align=True)
        row.operator("batoms.update", icon="FILE_REFRESH")
        row = layout.row()
        col = row.column()
        col.label(text="Dependencies:")
        #
        if not install_dependencies.has_module("ase"):
            op = layout.operator("batoms.install_package",
                                 icon='GREASEPENCIL', text="Install ase")
            op.package = "ase"
            op.modname = "ase"
        else:
            self.has_ase = True
            col.prop(self, "has_ase", text="ASE")
        #
        if not install_dependencies.has_module("skimage"):
            op = layout.operator("batoms.install_package",
                                 icon='GREASEPENCIL', text="Install scikit-image")
            op.package = "scikit-image"
            op.modname = "skimage"
        else:
            self.has_scikit_image = True
            col.prop(self, "has_scikit_image", text="scikit-image")
        #
        if not install_dependencies.has_module("spglib"):
            op = layout.operator("batoms.install_package",
                                 icon='GREASEPENCIL', text="Install spglib")
            op.package = "spglib"
            op.modname = "spglib"
        else:
            self.has_spglib = True
            col.prop(self, "has_spglib", text="spglib")
        # pymatgen
        if not install_dependencies.has_module("pymatgen"):
            op = layout.operator("batoms.install_package",
                                 icon='GREASEPENCIL', text="Install pymatgen")
            op.package = "pymatgen"
            op.modname = "pymatgen"
        else:
            self.has_pymatgen = True
            col.prop(self, "has_pymatgen", text="pymatgen")
        # openbabel
        if not install_dependencies.has_module("openbabel"):
            op = layout.operator("batoms.install_package",
                                 icon='GREASEPENCIL', text="Install openbabel")
            op.package = "openbabel"
            op.modname = "openbabel"
        else:
            self.has_openbabel = True
            col.prop(self, "has_openbabel", text="openbabel")

from . import install
    
classes = [BatomsAddonPreferences, 
install.BatomsUpdateButton]

def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)

def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)