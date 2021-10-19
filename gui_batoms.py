from ase.build import molecule, bulk
from ase import Atom, Atoms
import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (StringProperty,
                       BoolProperty,
                       BoolVectorProperty,
                       IntProperty,
                       IntVectorProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       PointerProperty,
                       )
from batoms.gui_io import import_batoms
from batoms.butils import get_selected_batoms
from batoms.batoms import Batoms

# The panel.
class Batoms_PT_prepare(Panel):
    bl_label       = "Batoms Tools"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    # bl_options     = {}
    bl_category = "Batoms"
    bl_idname = "BATOMS_PT_Tools"

    def draw(self, context):
        layout = self.layout
        bapanel = context.scene.bapanel

        box = layout.box()
        col = box.column()
        col.label(text="Model type")
        col = box.column()
        col.prop(bapanel, "model_type", expand  = True)
        box = layout.box()
        col = box.column()
        col.label(text="Polyhedra type")
        row = box.row()
        row.prop(bapanel, "polyhedra_type", expand  = True)
        row = box.row()
        row.prop(bapanel, "hide", expand  = True)
        row = box.row()
        row.prop(bapanel, "scale")

        box = layout.box()
        col = box.column(align=True)
        col.label(text="Add structure")
        col.prop(bapanel, "atoms_str")
        col.prop(bapanel, "atoms_name")
        # col = box.column()
        col.operator("batoms.add_molecule")
        col.operator("batoms.add_bulk")
        col.operator("batoms.add_atoms")

        box = layout.box()
        col = box.column(align=True)
        col.label(text="Export atoms")
        col.prop(bapanel, "filetype")

class BatomsProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    def Callback_model_type(self, context):
        bapanel = bpy.context.scene.bapanel
        model_type = list(bapanel.model_type)[0]
        modify_batoms_attr(self.selected_batoms, 'model_type', model_type)
    def Callback_polyhedra_type(self, context):
        bapanel = bpy.context.scene.bapanel
        polyhedra_type = list(bapanel.polyhedra_type)[0]
        modify_batoms_attr(self.selected_batoms, 'polyhedra_type', polyhedra_type)
    def Callback_modify_hide(self, context):
        bapanel = bpy.context.scene.bapanel
        modify_batoms_attr(self.selected_batoms, 'hide', bapanel.hide)
    def Callback_modify_scale(self, context):
        bapanel = bpy.context.scene.bapanel
        modify_batoms_attr(self.selected_batoms, 'scale', bapanel.scale)
    def Callback_export_atoms(self, context):
        bapanel = bpy.context.scene.bapanel
        export_atoms(self.selected_batoms, bapanel.filetype)
    
    model_type: EnumProperty(
        name="model_type",
        description="Structural models",
        items=(('0',"Space-filling", "Use ball and stick"),
               ('1',"Ball-and-stick", "Use ball"),
               ('2',"Polyhedral","Use polyhedral"),
               ('3',"Stick", "Use stick")),
        default={'0'}, 
        update=Callback_model_type,
        options={'ENUM_FLAG'},
        )
    polyhedra_type: EnumProperty(
        name="polyhedra_type",
        description="Polhhedra models",
        items=(('0',"0", "atoms, bonds and polyhedra"),
               ('1',"1", "atoms, polyhedra"),
               ('2',"2","central atoms, polyhedra"),
               ('3',"3", "polyhedra")),
        default={'0'}, 
        update=Callback_polyhedra_type,
        options={'ENUM_FLAG'},
        )
    hide: BoolProperty(name="hide", default=False, 
                update = Callback_modify_hide)
    scale: FloatProperty(
        name="scale", default=1.0,
        description = "scale", update = Callback_modify_scale)
    atoms_str: StringProperty(
        name = "Formula", default='H2O',
        description = "atoms_str")
    atoms_name: StringProperty(
        name = "Label", default='h2o',
        description = "Label")
    filetype: StringProperty(
        name = "File type", default='xyz',
        description = "save batoms to file", update = Callback_export_atoms)


def modify_batoms_attr(batoms_name_list, key, value):
    batoms_list = []
    for name in batoms_name_list:
        batoms = Batoms(label = name)
        setattr(batoms, key, value)
        batoms_list.append(batoms)
    for batoms in batoms_list:
        batoms.select = True
def export_atoms(batoms_name_list, filetype = 'xyz'):
    batoms_list = []
    for name in batoms_name_list:
        batoms = Batoms(label = name)
        batoms.write('%s.%s'%(batoms.label, filetype))
        batoms_list.append(batoms)
    for batoms in batoms_list:
        batoms.select = True


class AddMolecule(Operator):
    bl_idname = "batoms.add_molecule"
    bl_label = "Add molecule"
    bl_description = ("Add molecule")
    def execute(self, context):
        #print(molecule_str)
        bapanel = context.scene.bapanel
        atoms = molecule(bapanel.atoms_str)
        Batoms(label = bapanel.atoms_name, atoms = atoms)
        return {'FINISHED'}
class AddBulk(Operator):
    bl_idname = "batoms.add_bulk"
    bl_label = "Add bulk"
    bl_description = ("Add bulk")
    def execute(self, context):
        bapanel = context.scene.bapanel
        atoms = bulk(bapanel.atoms_str)
        import_batoms(atoms, name = bapanel.atoms_name)
        return {'FINISHED'}     
class AddAtoms(Operator):
    bl_idname = "batoms.add_atoms"
    bl_label = "Add atoms"
    bl_description = ("Add atoms")
    def execute(self, context):
        bapanel = context.scene.bapanel
        atoms = Atoms(bapanel.atoms_str)
        import_batoms(atoms, name = bapanel.atoms_name)
        return {'FINISHED'}
