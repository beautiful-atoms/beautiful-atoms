import bpy
from bpy.types import Panel
from bpy.props import (BoolProperty,
                       IntProperty,
                       FloatProperty,
                       EnumProperty,
                       )

from batoms.utils.butils import get_selected_vertices
from batoms.batoms import Batoms
from batoms.bond import Bond
from batoms.bond.slicebonds import SliceBonds
from batoms.gui.utils import get_attr, set_attr, get_enum_attr

# The panel.
class Bond_PT_prepare(Panel):
    bl_label = "Edit bonds"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Bond"
    bl_idname = "BATOMS_PT_Bond"
    bl_parent_id = 'VIEW3D_PT_Batoms_bond'

    # @classmethod
    # def poll(cls, context):
    #     obj = context.object
    #     if obj:
    #         return obj.batoms.type == 'BOND' and obj.mode != 'EDIT'
    #     else:
    #         return False

    def draw(self, context):
        layout = self.layout
        bond = context.scene.Bbond

        layout.label(text="Bond style")
        col = layout.column()
        col.prop(bond, "style", expand=True)
        layout.prop(bond, "order")
        # layout.prop(bond, "width")
        layout.prop(bond, "show")

# ---------------------------------------------------
def get_active_bond():
    context = bpy.context
    if context.object and context.object.batoms.type == 'BOND':
        indices = get_selected_vertices(context.object)
        if len(indices) > 0:
            batoms = Batoms(label=context.object.batoms.label)
            slicebonds = batoms.bond[indices]
            return slicebonds
        else:
            return None
    return None

def set_bond_attr_by_batoms(key):
    """
    """
    def setter(self, value):
        slicebonds = get_active_bond()
        if slicebonds is not None:
            # batoms = Batoms(label=slicebonds.label)
            # slicebonds = batoms.bond[slicebonds.indices]
            setattr(slicebonds, key, value)
            # bpy.ops.object.mode_set(mode="EDIT")
            bpy.context.view_layer.objects.active = slicebonds.obj

    return setter

def get_bond_attr(key):
    """
    """
    def getter(self):
        bond = get_active_bond()
        if bond is not None:
            # batoms = Batoms(label=bond.label)
            # bpy.context.view_layer.objects.active = batoms.bond.obj
            value = getattr(bond, key)[0]
            # bpy.ops.object.mode_set(mode="EDIT")
            return value
        else:
            return self.bl_rna.properties[key].default

    return getter


class BondProperties(bpy.types.PropertyGroup):


    style: EnumProperty(
        name="style",
        description="bond style",
        items=(("Unicolor cylinder", "Unicolor", "", 0),
               ("Bicolor cylinder", "Bicolor", "", 1),
               ("Dashed line", "Dashed", "", 2),
               ("Spring", "Spring", "", 3)
               ),
        default=1,
        get=get_enum_attr("style", get_active_bond),
        set=set_bond_attr_by_batoms("style"),
    )

    width: FloatProperty(
        name="bondwidth", default=0.1,
        min=0, soft_max=1,
        description="bondwidth",
        get=get_bond_attr("width"),
        set=set_bond_attr_by_batoms("width"),
    )

    order: IntProperty(name="Bond order", default=1,
                       min=1, max=3,
                       get=get_bond_attr("order"),
                       set=set_bond_attr_by_batoms("order"),
                       )
    show: BoolProperty(name="show", default=False,
                       get=get_bond_attr("show"),
                       set=set_bond_attr_by_batoms("show"),
                       )
