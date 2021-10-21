import bpy
from bpy.types import (Panel,
                       Operator,
                       AddonPreferences,
                       PropertyGroup,
                       )
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       PointerProperty,
                       )

from batoms.butils import get_selected_objects, get_selected_batoms
from batoms import Batoms

# The panel.
class Bond_PT_prepare(Panel):
    bl_label       = "Bond"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    # bl_options     = {'DEFAULT_CLOSED'}
    bl_category = "Bond"
    bl_idname = "BATOMS_PT_Bond"

  
    def draw(self, context):
        layout = self.layout
        bbpanel = context.scene.bbpanel

        box = layout.box()
        col = box.column()
        col.label(text="Bond style")
        col = box.column()
        col.prop(bbpanel, "bond_style", expand  = True)
        box = layout.box()
        col = box.column(align=True)
        row = box.row()
        row.prop(bbpanel, "min")
        row = box.row()
        row.prop(bbpanel, "max")
        row = box.row()
        row.prop(bbpanel, "bondwidth")

        box = layout.box()
        row = box.row()
        row.prop(bbpanel, "order")
        row = box.row()
        row.prop(bbpanel, "order_offset")

        box = layout.box()
        row = box.row()
        row.prop(bbpanel, "search")
        row = box.row()
        row.prop(bbpanel, "polyhedra")

        box = layout.box()
        col = box.column(align=True)
        row = box.row()
        row.prop(bbpanel, "bondcolor")

class BondProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    @property
    def selected_bond(self):
        return get_selected_objects('bbond')
    def Callback_bond_style(self, context):
        bbpanel = bpy.context.scene.bbpanel
        bond_style = list(bbpanel.bond_style)[0]
        modify_bond_attr(self.selected_batoms, self.selected_bond, 'style', bond_style)
    def Callback_modify_min(self, context):
        bbpanel = bpy.context.scene.bbpanel
        min = bbpanel.min
        modify_bond_attr(self.selected_batoms, self.selected_bond, 'min', min)
    def Callback_modify_max(self, context):
        bbpanel = bpy.context.scene.bbpanel
        max = bbpanel.max
        modify_bond_attr(self.selected_batoms, self.selected_bond, 'max', max)
    def Callback_modify_bondwidth(self, context):
        bbpanel = bpy.context.scene.bbpanel
        bondwidth = bbpanel.bondwidth
        modify_bond_attr(self.selected_batoms, self.selected_bond, 'width', bondwidth)
    def Callback_modify_search(self, context):
        bbpanel = bpy.context.scene.bbpanel
        search = bbpanel.search
        modify_bond_attr(self.selected_batoms, self.selected_bond, 'search', search)
    def Callback_modify_polyhedra(self, context):
        bbpanel = bpy.context.scene.bbpanel
        polyhedra = bbpanel.polyhedra
        modify_bond_attr(self.selected_batoms, self.selected_bond, 'polyhedra', polyhedra)
    def Callback_modify_bondcolor(self, context):
        bbpanel = bpy.context.scene.bbpanel
        bondcolor = bbpanel.bondcolor
        modify_bond_attr(self.selected_batoms, self.selected_bond, 'color', bondcolor)
    def Callback_modify_order(self, context):
        bbpanel = bpy.context.scene.bbpanel
        order = bbpanel.order
        modify_bond_attr(self.selected_batoms, self.selected_bond, 'order', order)
    def Callback_modify_order_offset(self, context):
        bbpanel = bpy.context.scene.bbpanel
        order_offset = bbpanel.order_offset
        modify_bond_attr(self.selected_batoms, self.selected_bond, 'order_offset', order_offset)

    bond_style: EnumProperty(
        name="style",
        description="bond style",
        items=(('0',"Unicolor cylinder", ""),
               ('1',"Bicolor cylinder", ""),
               ('2',"Dashed line", ""),
               ('3',"Dotted line", "")),
        default={'1'},
        update=Callback_bond_style,
        options={'ENUM_FLAG'},
        )
    bondwidth: FloatProperty(
        name="bondwidth", default=0.1,
        description = "bondwidth", update = Callback_modify_bondwidth)
    min: FloatProperty(
        name="Length min", default=0,
        description = "min", update = Callback_modify_min)
    max: FloatProperty(
        name="Length max", default=2.0,
        description = "max", update = Callback_modify_max)
    search: IntProperty(name="Search mode", default=0, 
                update = Callback_modify_search)
    polyhedra: BoolProperty(name="polyhedra", default=False, 
                update = Callback_modify_polyhedra)
    bondcolor: FloatVectorProperty(
        name="bondcolor", 
        subtype='COLOR',
        default=(0.1, 0.8, 0.4 ,1.0),
        size =4,
        description="color picker",
        update = Callback_modify_bondcolor)
    order: IntProperty(name="Bond order", default=1, 
                update = Callback_modify_order)
    order_offset: FloatProperty(name="order_offset", default=0.15, 
                update = Callback_modify_order_offset)

def modify_bond_attr(selected_batoms, selected_bond, key, value):
    selected_bond_new = []
    for batoms_name in selected_batoms:
        batoms = Batoms(label = batoms_name)
        for bond_name in selected_bond:
            bond = bpy.data.objects[bond_name]
            if bond.bbond.label == batoms_name:
                if key != 'color':
                    setattr(batoms.bondsetting['%s-%s'%(bond.bbond.species1, 
                        bond.bbond.species2)], key, value)
                else:
                    index = [bond.bbond.species1, bond.bbond.species2].index(bond.bbond.species) + 1
                    setattr(batoms.bondsetting['%s-%s'%(bond.bbond.species1, 
                            bond.bbond.species2)], 'color%s'%index, value)
                if batoms.bondsetting['%s-%s'%(bond.bbond.species1, bond.bbond.species2)].style == '0':
                    selected_bond_new.append('bond_%s_%s_%s'%(bond.bbond.label, 
                        bond.bbond.species1, bond.bbond.species2))
                else:
                    selected_bond_new.append('bond_%s_%s_%s_%s'%(bond.bbond.label, 
                        bond.bbond.species1, bond.bbond.species2, bond.bbond.species))
        batoms.draw_bonds()
        if batoms.model_type == 2:
            batoms.draw_polyhedras()
    for name in selected_bond_new:
        obj = bpy.data.objects.get(name)
        obj.select_set(True)        