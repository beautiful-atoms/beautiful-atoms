import bpy
from bpy.types import Panel
from bpy.props import (BoolProperty,
                       IntProperty,
                       FloatProperty,
                       EnumProperty,
                       )

from batoms.utils.butils import get_selected_vertices
from batoms.batoms import Batoms
from batoms.bond.bond import Bond
from batoms.gui.gui_batoms import get_attr, set_attr, get_enum_attr

# The panel.
class Bond_PT_prepare(Panel):
    bl_label = "Bond"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Bond"
    bl_idname = "BATOMS_PT_Bond"

    # @classmethod
    # def poll(cls, context):
    #     obj = context.object
    #     if obj:
    #         return obj.batoms.type == 'BOND' and obj.mode != 'EDIT'
    #     else:
    #         return False
    
    def draw(self, context):
        layout = self.layout
        bond = context.scene.batoms.bond

        layout.operator("bond.bond_order_auto_set", icon='MODIFIER_ON', text="Auto Set Bond Order")
        layout.separator()
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
        v = get_selected_vertices(context.object)
        if len(v) > 0:
            bond = Bond(label=context.object.batoms.label, index=v[0])
            return bond
        else:
            return None
    return None

def set_bond_attr_by_batoms(key):
    """
    """
    def setter(self, value):
        bond = get_active_bond()
        if bond is not None:
            batoms = Batoms(label=bond.label)
            setattr(batoms.bonds[bond.index], key, value)
            # bpy.ops.object.mode_set(mode="EDIT")
            bpy.context.view_layer.objects.active = batoms.bonds.obj
    
    return setter

def get_bond_attr(key):
    """
    """
    def getter(self):
        bond = get_active_bond()
        if bond is not None:
            batoms = Batoms(label=bond.label)
            # bpy.context.view_layer.objects.active = batoms.bonds.obj
            value = getattr(batoms.bonds[bond.index], key)
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
                       get=get_attr("order", get_active_bond),
                       set=set_bond_attr_by_batoms("order"),
                       )
    show: BoolProperty(name="show", default=False,
                       get=get_attr("show", get_active_bond),
                       set=set_bond_attr_by_batoms("show"),    
                       )
