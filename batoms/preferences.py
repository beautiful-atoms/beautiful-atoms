"""

1) set custom file folder.
2) use default startup of batoms.
3) use default preferences of batoms.
#TODO  batoms_setting_path_update
"""

import bpy
from bpy.types import AddonPreferences
from bpy.props import (
    BoolProperty,
    StringProperty,
)
from batoms.install.pip_dependencies import has_module
from batoms.install import update

dependencies = {"ase": "ase",
                "scikit-image": "skimage",
                "spglib": "spglib",
                "pymatgen": "pymatgen",
                "openbabel": "openbabel",
                }


DEFAULT_GITHUB_ACCOUNT = "superstar54"
DEFAULT_REPO_NAME = "beautiful-atoms"
DEFAULT_PLUGIN_NAME = "batoms"


class BatomsDefaultPreference(bpy.types.Operator):
    """Update Batoms"""
    bl_idname = "batoms.use_batoms_preference"
    bl_label = "Use defatul preference of Batoms"
    bl_description = "Use startup file of Batoms"

    def execute(self, context):
        bpy.context.preferences.view.use_translate_new_dataname = False
        bpy.context.preferences.inputs.use_rotate_around_active = True
        bpy.context.preferences.inputs.use_zoom_to_mouse = True
        # For laptop
        bpy.context.preferences.inputs.use_emulate_numpad = True
        # For laptop without mouse
        bpy.context.preferences.inputs.use_mouse_emulate_3_button = True
        # bpy.context.window.workspace = bpy.data.workspaces['UV Editing']
        self.report({"INFO"}, "Set default preferences successfully!")
        return {'FINISHED'}

class BatomsDefaultStartup(bpy.types.Operator):
    """Update Batoms"""
    bl_idname = "batoms.use_batoms_startup"
    bl_label = "Use startup file of Batoms"
    bl_description = "Use defatul startup of Batoms"

    def execute(self, context):
        import sys, os
        import pathlib
        addon_dir = pathlib.Path(__file__).parent.resolve()
        blend_dir = os.path.join(addon_dir, "data/startup.blend")
        bpy.ops.wm.open_mainfile(filepath=blend_dir, load_ui=True, use_scripts=True) 
        bpy.ops.wm.save_homefile()
        self.report({"INFO"}, "Load default startup successfully!")
        # todo open preference again.
        # bpy.ops.screen.userpref_show('INVOKE_DEFAULT')
        # bpy.ops.preferences.addon_show(module="batoms")
        return {'FINISHED'}

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
        #
        layout.operator("batoms.use_batoms_startup", icon="FILE_REFRESH")
        layout.operator("batoms.use_batoms_preference", icon="FILE_REFRESH")
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

classes = [BatomsDefaultPreference,
            BatomsDefaultStartup,
            BatomsAddonPreferences,
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

