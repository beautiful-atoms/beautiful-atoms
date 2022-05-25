"""

0) update batoms
1) set custom file folder.
2) set custom windows

"""

import bpy
from bpy.types import AddonPreferences
from bpy.props import (
    BoolProperty,
    StringProperty,
)
from batoms.install.pip_dependencies import has_module
import pathlib
from batoms.install import update

dependencies = {"ase": "ase",
                "scikit-image": "skimage",
                "spglib": "spglib",
                "pymatgen": "pymatgen",
                "openbabel": "openbabel",
                }


# bpy.context.window.workspace = bpy.data.workspaces['UV Editing']

class BatomsDefaultStartup(bpy.types.Operator):
    """Update Batoms"""
    bl_idname = "batoms.use_batoms_startup"
    bl_label = "Use startup file of Batoms"
    bl_description = "Use startup file of Batoms"

    def execute(self, context):
        if not has_git():
            self.report({"ERROR"}, "Please install Git first.")
            print("Please install Git first.")
            return {"CANCELLED"}
        gitclone()
        return {"FINISHED"}

class BatomsAddonPreferences(AddonPreferences):
    bl_idname = __package__

    def batoms_setting_path_update(self, context):
        import os
        import subprocess
        if os.name == 'posix':  # Linux
            cmds = ["export", "BATOMS_SETTING_PATH={}".format(
                self.batoms_setting_path)]
        if os.name == 'nt':  # Windows
            cmds = ["setx", "BATOMS_SETTING_PATH {}".format(
                self.batoms_setting_path)]
        print(update.run(cmds))

    ase: BoolProperty(
        name="ASE installed",
        description="ASE package installed",
        default=False,
    )

    skimage: BoolProperty(
        name="scikit-image installed",
        description="scikit-image package installed",
        default=False,
    )

    spglib: BoolProperty(
        name="spglib installed",
        description="spglib package installed",
        default=False,
    )

    pymatgen: BoolProperty(
        name="pymatgen installed",
        description="pymatgen package installed",
        default=False,
    )

    openbabel: BoolProperty(
        name="Openbabel installed",
        description="openbabel package installed",
        default=False,
    )

    batoms_setting_path: StringProperty(name="Custom Setting Path", description="Custom Setting Path",
                                            default="", subtype="FILE_PATH", update=batoms_setting_path_update)


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
        layout.separator()
        row = layout.row()
        col = row.column()
        col.label(text="Dependencies:")
        #
        for package, modname in dependencies.items():
            if not has_module(modname):
                op = layout.operator("batoms.pip_install_package",
                                     icon='GREASEPENCIL', text="Install {}".format(package))
                op.package = package
                op.modname = modname
            else:
                setattr(self, modname, True)
                col.prop(self, modname, text=package)
        # custom folder
        layout.separator()
        box = layout.box().column()
        box.label(text="Custom Settings")
        box.prop(self, "batoms_setting_path")

classes = [BatomsAddonPreferences,
           update.BatomsUpdateButton,
           ]


def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)


def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
