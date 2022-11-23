"""
#TODO wrap is 3d vector, do we need support all?
"""
import bpy
from bpy.types import Panel
from bpy.props import (
    BoolProperty,
    FloatProperty,
    EnumProperty,
    StringProperty,
)
from batoms import Batoms
from batoms.gui.utils import (get_active_bpy_data,
        get_attr, get_enum_attr, set_attr, set_enum_attr)

model_style_items = [("Space-filling", "Space-filling", "", 0),
                     ("Ball-and-stick", "Ball-and-Stick", "", 1),
                     ("Polyhedral", "Polyhedral", "", 2),
                     ("Stick", "Stick", "", 3),
                     ]


class Batoms_PT_prepare(Panel):
    bl_label = "Batoms"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "BATOMS_PT_Tools"

    def draw(self, context):
        name = 'None'
        if context.object:
            if context.object.batoms.type != 'OTHER':
                name = context.object.batoms.label
        layout = self.layout
        layout.label(text="Active: " + name)
        batoms = context.scene.batoms.batoms

        layout.label(text="Model style")
        col = layout.column()
        col.prop(batoms, "model_style", expand=True)
        # layout.label(text="Add label")
        layout.label(text="Radius style")
        layout.prop(batoms, "radius_style", expand=True)
        layout.label(text="Color style")
        layout.prop(batoms, "color_style", expand=True)
        layout.label(text="Polyhedra style")
        layout.prop(batoms, "polyhedra_style", expand=True)

        split = layout.split()
        col = split.column()
        col.prop(batoms, "show", expand=True)
        col = split.column()
        col.prop(batoms, "wrap", expand=True)
        layout.prop(batoms, "crystal_view", expand=True)
        layout.prop(batoms, "show_label", expand=True, text="label")
        layout.prop(batoms, "scale")

        # layout.operator("batoms.replace")





# ---------------------------------------------------
def get_active_batoms():
    context = bpy.context
    if context.object and context.object.batoms.type != 'OTHER':
        mode = context.object.mode
        batoms = Batoms(label=context.object.batoms.label)
        bpy.context.view_layer.objects.active = batoms.obj
        bpy.ops.object.mode_set(mode=mode)
        return batoms
    return None


def set_batoms_attr(key, value):
    """
    """
    batoms = get_active_batoms()
    if batoms is not None:
        setattr(batoms, key, value)
        bpy.context.view_layer.objects.active = batoms.obj


def modify_batoms_attr(context, key, value):
    """
    """
    if context.object and context.object.batoms.type != 'OTHER':
        batoms = Batoms(label=context.object.batoms.label)
        setattr(batoms, key, value)
        # batoms.obj.select_set(True)
        bpy.context.view_layer.objects.active = batoms.obj


def get_wrap(self):
    batoms = get_active_bpy_data('batoms')()
    if batoms is not None:
        return batoms.wrap[0]
    else:
        return 0


class BatomsProperties(bpy.types.PropertyGroup):
    model_style: EnumProperty(
        name="model_style",
        description="Structural models",
        items=model_style_items,
        get=get_enum_attr("model_style", get_active_bpy_data('batoms')),
        set=set_enum_attr("model_style", set_batoms_attr),
        default=0,
    )

    show_label: StringProperty(
        name="label",
        description="Show label: None, Index, Species or Charge and so on",
        get=get_attr("show_label", get_active_bpy_data('batoms')),
        set=set_attr("show_label", set_batoms_attr),
        default="",
    )

    radius_style: EnumProperty(
        name="radius_style",
        description="Structural models",
        items=(("Covalent", "covalent", "", 0),
               ("VDW", "van der Waals", "", 1),
               ("Ionic", "ionic", "", 2)),
        get=get_enum_attr("radius_style", get_active_bpy_data('batoms')),
        set=set_enum_attr("radius_style", set_batoms_attr),
        default=0,
    )

    color_style: EnumProperty(
        name="color_style",
        description="Color",
        items=(("JMOL", "JMOL", "", 0),
               ("VESTA", "VESTA", "", 1),
               ("CPK", "CPK", "", 2)),
        get=get_enum_attr("color_style", get_active_bpy_data('batoms')),
        set=set_enum_attr("color_style", set_batoms_attr),
        default=0,
    )

    polyhedra_style: EnumProperty(
        name="polyhedra_style",
        description="Polhhedra models",
        items=(("0", "atoms, bonds and polyhedra", "", 0),
               ("1", "atoms, polyhedra", "", 1),
               ("2", "central atoms, polyhedra", "", 2),
               ("3", "polyhedra", "", 3)),
        get=get_enum_attr("polyhedra_style", get_active_bpy_data('batoms')),
        set=set_enum_attr("polyhedra_style", set_batoms_attr),
        default=0,
    )

    show: BoolProperty(name="show",
                       default=False,
                       description="show all object for view and rendering",
                       get=get_attr("show", get_active_bpy_data('batoms')),
                       set=set_attr("show", set_batoms_attr)
                       )

    wrap: BoolProperty(name="wrap",
                       default=False,
                       description="wrap all atoms into cell",
                       get=get_wrap,
                       set=set_attr("wrap", set_batoms_attr)
                       )

    crystal_view: BoolProperty(name="crystal_view",
                       default=False,
                       description="crystal_view all object for view and rendering",
                       get=get_attr("crystal_view", get_active_bpy_data('batoms')),
                       set=set_attr("crystal_view", set_batoms_attr)
                       )

    scale: FloatProperty(
        name="scale", default=1.0,
        min=0.0, soft_max=2.0,
        description="scale",
        get=get_attr("scale", get_active_bpy_data('batoms')),
        set=set_attr("scale", set_batoms_attr)
    )
