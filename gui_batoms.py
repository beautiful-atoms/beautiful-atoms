from ase.build import molecule, bulk
from ase import Atom, Atoms
import bpy
import bmesh
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
from batoms.butils import read_batoms_select
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
        blpanel = context.scene.blpanel

        box = layout.box()
        col = box.column()
        col.label(text="Model")
        row = box.row()
        row.prop(blpanel, "model_type", expand  = True)
        row = box.row()
        row.prop(blpanel, "scale")

        box = layout.box()
        col = box.column(align=True)
        col.label(text="Add structure")
        col.prop(blpanel, "atoms_str")
        col.prop(blpanel, "atoms_name")
        # col = box.column()
        col.operator("batoms.add_molecule")
        col.operator("batoms.add_bulk")
        col.operator("batoms.add_atoms")

        box = layout.box()
        col = box.column(align=True)
        col.label(text="Export atoms")
        col.prop(blpanel, "filetype")



class BatomsProperties(bpy.types.PropertyGroup):
    @property
    def batoms_list(self):
        return self.get_batoms_list()
    def get_batoms_list(self):
        batoms_list = read_batoms_select()
        return batoms_list
    def Callback_model_type(self, context):
        blpanel = bpy.context.scene.blpanel
        print('Callback_model_type')
        model_type = list(blpanel.model_type)[0]
        modify_model_type(self.batoms_list, model_type)
    def Callback_modify_scale(self, context):
        blpanel = bpy.context.scene.blpanel
        print('Callback_modify_scale')
        modify_scale(self.batoms_list, blpanel.scale)
    def Callback_export_atoms(self, context):
        blpanel = bpy.context.scene.blpanel
        print('Callback_export_atoms')
        export_atoms(self.batoms_list, blpanel.filetype)
    
    model_type: EnumProperty(
        name="Model",
        description="Structural models",
        items=(('0',"Space-filling", "Use ball and stick"),
               ('1',"Ball-and-stick", "Use ball"),
               ('2',"Polyhedral","Use polyhedral"),
               ('3',"Stick", "Use stick")),
        default={'0'}, 
        update=Callback_model_type,
        options={'ENUM_FLAG'},
        )
    scale: FloatProperty(
        name="scale", default=1.0,
        description = "scale", update = Callback_modify_scale)
    atoms_str: StringProperty(
        name = "Formula", default='H2O',
        description = "atoms_str")
    atoms_name: StringProperty(
        name = "Name", default='h2o',
        description = "name")
    filetype: StringProperty(
        name = "File type", default='xyz',
        description = "save batoms to file", update = Callback_export_atoms)
    

class AddMolecule(Operator):
    bl_idname = "batoms.add_molecule"
    bl_label = "Add molecule"
    bl_description = ("Add molecule")
    def execute(self, context):
        #print(molecule_str)
        blpanel = context.scene.blpanel
        atoms = molecule(blpanel.atoms_str)
        Batoms(label = blpanel.atoms_name, atoms = atoms)
        return {'FINISHED'}
class AddBulk(Operator):
    bl_idname = "batoms.add_bulk"
    bl_label = "Add bulk"
    bl_description = ("Add bulk")
    def execute(self, context):
        blpanel = context.scene.blpanel
        atoms = bulk(blpanel.atoms_str)
        import_batoms(atoms, name = blpanel.atoms_name)
        return {'FINISHED'}     
class AddAtoms(Operator):
    bl_idname = "batoms.add_atoms"
    bl_label = "Add atoms"
    bl_description = ("Add atoms")
    def execute(self, context):
        blpanel = context.scene.blpanel
        atoms = Atoms(blpanel.atoms_str)
        import_batoms(atoms, name = blpanel.atoms_name)
        return {'FINISHED'}     

def modify_model_type(batoms_name_list, model_type):
    for name in batoms_name_list:
        batoms = Batoms(label = name)
        batoms.model_type = model_type
def modify_scale(batoms_name_list, scale):
    for name in batoms_name_list:
        batoms = Batoms(label = name)
        for species, ba in batoms.batoms.items():
            ba.scale = scale
def export_atoms(batoms_name_list, filetype = 'xyz'):
    for name in batoms_name_list:
        batoms = Batoms(label = name)
        batoms.write('%s.%s'%(batoms.label, filetype))