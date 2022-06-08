import bpy
from bpy.types import Panel
from bpy.props import (
    BoolProperty,
    FloatProperty,
    EnumProperty,
    StringProperty,
)
from batoms import Batoms


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
        layout.operator("batoms.import")
        layout.operator("batoms.export")

        bapanel = context.scene.bapanel

        layout.label(text="Model style")
        col = layout.column()
        col.prop(bapanel, "model_style", expand=True)
        # layout.label(text="Add label")
        layout.label(text="Radius style")
        layout.prop(bapanel, "radius_style", expand=True)
        layout.label(text="Color style")
        layout.prop(bapanel, "color_style", expand=True)
        layout.label(text="Polyhedra style")
        layout.prop(bapanel, "polyhedra_style", expand=True)

        layout.prop(bapanel, "show", expand=True)
        layout.prop(bapanel, "wrap", expand=True)
        layout.prop(bapanel, "show_label", expand=True, text="label")
        layout.prop(bapanel, "scale")

        layout.operator("batoms.replace")


def get_model_style(self):
    batoms = get_active_collection()
    if batoms is not None:
        return int(batoms.model_style)
    else:
        return 0

def set_model_style(self, value):
    items = self.bl_rna.properties["model_style"].enum_items
    item = items[value]
    model_style = item.identifier
    # print(model_style)
    self["model_style"] = model_style
    set_batoms_attr('model_style', value)


def get_radius_style(self):
    batoms = get_active_collection()
    if batoms is not None:
        return int(batoms.radius_style)
    else:
        return 0

def set_radius_style(self, value):
    items = self.bl_rna.properties["radius_style"].enum_items
    item = items[value]
    radius_style = item.identifier
    # print(radius_style)
    self["radius_style"] = radius_style
    set_batoms_attr('radius_style', value)


def get_color_style(self):
    batoms = get_active_collection()
    if batoms is not None:
        return int(batoms.color_style)
    else:
        return 0

def set_color_style(self, value):
    items = self.bl_rna.properties["color_style"].enum_items
    item = items[value]
    color_style = item.identifier
    # print(color_style)
    self["color_style"] = color_style
    set_batoms_attr('color_style', value)


def get_polyhedra_style(self):
    batoms = get_active_collection()
    if batoms is not None:
        return int(batoms.polyhedra_style)
    else:
        return 0

def set_polyhedra_style(self, value):
    items = self.bl_rna.properties["polyhedra_style"].enum_items
    item = items[value]
    polyhedra_style = item.identifier
    # print(polyhedra_style)
    self["polyhedra_style"] = polyhedra_style
    set_batoms_attr('polyhedra_style', value)

def get_show_label(self):
    batoms = get_active_collection()
    if batoms is not None:
        return batoms.show_label
    else:
        return ""

def set_show_label(self, value):
    self["show_label"] = value
    set_batoms_attr('show_label', value)

def get_show(self):
    batoms = get_active_collection()
    if batoms is not None:
        return batoms.show
    else:
        return False

def set_show(self, value):
    self["show"] = value
    set_batoms_attr('show', value)

def get_wrap(self):
    batoms = get_active_collection()
    if batoms is not None:
        return batoms.wrap[0]
    else:
        return 0

def set_wrap(self, value):
    self["wrap"] = value
    set_batoms_attr('wrap', value)

def get_scale(self):
    batoms = get_active_collection()
    if batoms is not None:
        return batoms.scale
    else:
        return 0

def set_scale(self, value):
    self["scale"] = value
    set_batoms_attr('scale', value)

class BatomsProperties(bpy.types.PropertyGroup):
    model_style: EnumProperty(
        name="model_style",
        description="Structural models",
        items=model_style_items,
        get=get_model_style,
        set=set_model_style,
        default=0,
    )

    show_label: StringProperty(
        name="label",
        description="Show label: None, Index, Species or Charge and so on",
        get=get_show_label,
        set=set_show_label,
        default="",
    )

    radius_style: EnumProperty(
        name="radius_style",
        description="Structural models",
        items=(("Covalent", "covalent", "", 0),
               ("VDW", "van der Waals", "", 1),
               ("Ionic", "ionic", "", 2)),
        get=get_radius_style,
        set=set_radius_style,
        default=0,
    )

    color_style: EnumProperty(
        name="color_style",
        description="Color",
        items=(("JMOL", "JMOL", "", 0),
               ("VESTA", "VESTA", "", 1),
               ("CPK", "CPK", "", 2)),
        get=get_color_style,
        set=set_color_style,
        default=0,
    )

    polyhedra_style: EnumProperty(
        name="polyhedra_style",
        description="Polhhedra models",
        items=(("0", "atoms, bonds and polyhedra", "", 0),
               ("1", "atoms, polyhedra", "", 1),
               ("2", "central atoms, polyhedra", "", 2),
               ("3", "polyhedra", "", 3)),
        get=get_polyhedra_style,
        set=set_polyhedra_style,
        default=0,
    )

    show: BoolProperty(name="show",
                       default=False,
                       description="show all object for view and rendering",
                       get=get_show,
                       set=set_show)

    wrap: BoolProperty(name="wrap",
                       default=False,
                       description="wrap all atoms into cell",
                       get=get_wrap,
                       set=set_wrap,)

    scale: FloatProperty(
        name="scale", default=1.0,
        min=0.0, soft_max=2.0,
        description="scale",
        get=get_scale,
        set=set_scale,)


def get_active_collection():
    """Get the collection of the active Batoms

    When get the attribute of Batoms object, 
    if the attribute if saved in the Batoms.coll.batoms,
    we only need to read data form the colleciton, 
    it is faster than get data from the Batoms itself.

    Returns:
        bpy.type.collection: _description_
    """
    context = bpy.context
    if context.object and context.object.batoms.type != 'OTHER':
        return bpy.data.collections[context.object.batoms.label].batoms
    return None

def get_active_batoms():
    context = bpy.context
    if context.object and context.object.batoms.type != 'OTHER':
        mode = context.object.mode
        batoms = Batoms(label=context.object.batoms.label)
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
