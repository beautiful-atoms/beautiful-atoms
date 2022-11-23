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
    EnumProperty,
)
from batoms.install.pip_dependencies import has_module
from batoms.install import update
from batoms.logger import update_logging_level
import logging

# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)


# Enum property.
# Note: the getter/setter callback must use integer identifiers!
logging_level_items = [
    ("DEBUG", "DEBUG", "", 0),
    ("INFO", "INFO", "", 1),
    ("WARNING", "WARNING", "", 2),
    ("ERROR", "ERROR", "", 3),
    ("CRITICAL", "CRITICAL", "", 4),
]

dependencies = {"ase": "ase",
                "scikit-image": "skimage",
                "spglib": "spglib",
                "pymatgen": "pymatgen",
                "openbabel": "openbabel",
                }


DEFAULT_GITHUB_ACCOUNT = "beautiful-atoms"
DEFAULT_REPO_NAME = "beautiful-atoms"
DEFAULT_PLUGIN_NAME = "batoms"

def get_plugin(key):
    """Helper function to get plugin
    Args:
        key (_type_): _description_
    """
    def getter(self):
        # self.get: returns the value of the custom property assigned to
        # key or default when not found
        value = self.bl_rna.properties[key].default
        return self.get(key, value)
    return getter

def set_plugin(key):
    """Helper function to set plugin

    Args:
        key (_type_): _description_
        value (_type_): _description_
    """
    def setter(self, value):
        import importlib
        self[key] = value
        plugin = importlib.import_module("batoms.plugins.{}".format(key))
        if value:
            plugin.register_class()
            logger.info("Enable {} plugin.".format(key))
        else:
            plugin.unregister_class()
            logger.info("Disnable {} plugin.".format(key))
    return setter

class BatomsDefaultPreference(bpy.types.Operator):
    """Update Batoms"""
    bl_idname = "batoms.use_batoms_preference"
    bl_label = "Use defatul preference of Batoms"
    bl_description = "Use startup file of Batoms"

    def execute(self, context):
        import pathlib
        import os
        batoms_asset_dir = os.path.join(pathlib.Path(__file__).parent.resolve(), 'asset')
        batoms_asset_dir = os.path.join(batoms_asset_dir, 'libraries')

        bpy.context.preferences.view.use_translate_new_dataname = False
        bpy.context.preferences.inputs.use_rotate_around_active = True
        bpy.context.preferences.inputs.use_zoom_to_mouse = True
        # For laptop
        bpy.context.preferences.inputs.use_emulate_numpad = True
        # For laptop without mouse
        bpy.context.preferences.inputs.use_mouse_emulate_3_button = True
        # bpy.context.window.workspace = bpy.data.workspaces['UV Editing']
        # theme
        bpy.context.preferences.themes[0].view_3d.space.gradients.background_type = "LINEAR"
        bpy.context.preferences.themes[0].view_3d.space.gradients.high_gradient = (0.9, 0.9, 0.9)
        bpy.context.preferences.themes[0].view_3d.space.gradients.gradient = (0.5, 0.5, 0.5)
        bpy.ops.wm.save_userpref()
        # logger
        bpy.context.preferences.addons['batoms'].preferences.logging_level = "WARNING"
        # asset_libraries
        if "Batoms" not in bpy.context.preferences.filepaths.asset_libraries.keys():
            bpy.ops.preferences.asset_library_add(directory=batoms_asset_dir)
            bpy.context.preferences.filepaths.asset_libraries[-1].name = "Batoms"
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
        ###################################################
        # Add additional settings to the startup file here
        ###################################################
        bpy.context.scene.unit_settings.system = "NONE"
        # viewport overlayrs
        for area in bpy.context.screen.areas:
            if area.type == 'VIEW_3D':
                for space in area.spaces:
                    if space.type == 'VIEW_3D':
                        space.overlay.show_extras = False
                        space.overlay.show_relationship_lines = False
                        break
        ###################################################
        bpy.ops.wm.save_homefile()
        self.report({"INFO"}, "Load default startup successfully!")
        # todo open preference again.
        # bpy.ops.screen.userpref_show('INVOKE_DEFAULT')
        # bpy.ops.preferences.addon_show(module="batoms")
        return {'FINISHED'}



class BatomsAddonPreferences(AddonPreferences):
    bl_idname = __package__

    def get_logging_level(self):
        items = self.bl_rna.properties["logging_level"].enum_items
        # self.get: returns the value of the custom property assigned to
        # key or default when not found
        return items[self.get("logging_level", 2)].value

    def set_logging_level(self, value):
        items = self.bl_rna.properties["logging_level"].enum_items
        item = items[value]
        level = item.identifier
        # we need to update both the preference and the logger
        self["logging_level"] = level
        # Set the logging level for all child loggers of "batoms"
        update_logging_level()
        # Note the following logging info might not emit
        # if global level is higher than INFO
        logger.info("Set logging level to: {}".format(level))



    def batoms_setting_path_update(self, context):
        import os
        import subprocess
        if os.name == 'posix':  # Linux
            cmds = ["export", "BATOMS_SETTING_PATH={}".format(
                self.batoms_setting_path)]
        if os.name == 'nt':  # Windows
            cmds = ["setx", "BATOMS_SETTING_PATH {}".format(
                self.batoms_setting_path)]
        logger.debug(update.subprocess_run(cmds))

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

    logging_level: EnumProperty(
        name="Logging Level",
        items=logging_level_items,
        get=get_logging_level,
        set=set_logging_level,
        default=2,
        )

    isosurface: BoolProperty(
        name="isosurface",
        description="Enable isosurface plugin",
        get=get_plugin("isosurface"),
        set=set_plugin("isosurface"),
        default=True,
    )

    molecular_surface: BoolProperty(
        name="molecular_surface",
        description="Enable molecular_surface plugin",
        get=get_plugin("molecular_surface"),
        set=set_plugin("molecular_surface"),
        default=True,
    )

    real_interaction: BoolProperty(
        name="real_interaction",
        description="Enable real_interaction plugin",
        get=get_plugin("real_interaction"),
        set=set_plugin("real_interaction"),
        default=False,
    )

    magres: BoolProperty(
        name="magres",
        description="Enable magres plugin",
        get=get_plugin("magres"),
        set=set_plugin("magres"),
        default=True,
    )

    highlight: BoolProperty(
        name="highlight",
        description="Enable highlight plugin",
        get=get_plugin("highlight"),
        set=set_plugin("highlight"),
        default=True,
    )

    cavity: BoolProperty(
        name="cavity",
        description="Enable cavity plugin",
        get=get_plugin("cavity"),
        set=set_plugin("cavity"),
        default=True,
    )

    crystal_shape: BoolProperty(
        name="crystal_shape",
        description="Enable crystal_shape plugin",
        get=get_plugin("crystal_shape"),
        set=set_plugin("crystal_shape"),
        default=True,
    )

    lattice_plane: BoolProperty(
        name="lattice_plane",
        description="Enable lattice_plane plugin",
        get=get_plugin("lattice_plane"),
        set=set_plugin("lattice_plane"),
        default=True,
    )

    template: BoolProperty(
        name="template",
        description="Enable template plugin",
        get=get_plugin("template"),
        set=set_plugin("template"),
        default=True,
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
        layout.label(text="Use default setting.")
        box = layout.box().column()
        row = box.row(align=True)
        row.operator("batoms.use_batoms_startup", text = "Use startup", icon="FILE_REFRESH")
        row.operator("batoms.use_batoms_preference", text = "Use Preferences", icon="FILE_REFRESH")
        layout.separator()

        box = layout.box().column()
        box.label(text="Dependencies:")
        #
        for package, modname in dependencies.items():
            if not has_module(modname):
                op = box.operator("batoms.pip_install_package",
                                     icon='IMPORT', text="Install {}".format(package))
                op.package = package
                op.modname = modname
            else:
                setattr(self, modname, True)
                box.prop(self, modname, text=package)
        #
        layout.separator()
        split = layout.split()
        col = split.column()
        col.label(text="Custom Plugins")
        col.prop(self, "highlight")
        col.prop(self, "isosurface")
        col.prop(self, "molecular_surface")
        col.prop(self, "crystal_shape")
        col.prop(self, "lattice_plane")
        col.prop(self, "cavity")
        col.prop(self, "magres")
        col.prop(self, "real_interaction")
        col.prop(self, "template")
        col = split.column()
        # custom folder
        layout.separator()
        box = layout.box().column()
        box.label(text="Custom Settings")
        box.prop(self, "logging_level")
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
