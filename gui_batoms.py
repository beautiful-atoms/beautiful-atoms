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
from batoms.butils import read_batoms_collection_list
from batoms.batoms import Batoms

# The panel.
class Batoms_PT_prepare(Panel):
    bl_label       = "Batoms Tools"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    # bl_options     = {}
    bl_category = "Batoms"
    bl_idname = "BATOMS_PT_Tools"

    filename: StringProperty(
        name = "Filename", default='batoms-output.xyz',
        description = "Export atoms to file.")

    def draw(self, context):
        layout = self.layout
        blpanel = context.scene.blpanel

        box = layout.box()
        col = box.row()
        col.prop(blpanel, "collection_list")
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
        col.label(text="Render atoms")
        col.prop(blpanel, "output_image")



class BatomsProperties(bpy.types.PropertyGroup):
    def Callback_model_type(self, context):
        blpanel = bpy.context.scene.blpanel
        print('Callback_model_type')
        model_type = list(blpanel.model_type)[0]
        modify_model_type(blpanel.collection_list, model_type)
    
    def Callback_collection_list(self, context):
        print('Callback_collection_list')
        items = read_batoms_collection_list()
        items = [(item, item, "") for item in items]
        items = tuple(items)
        return items
    def Callback_modify_scale(self, context):
        blpanel = bpy.context.scene.blpanel
        print('Callback_modify_scale')
        modify_scale(blpanel.collection_list, blpanel.scale)
    def Callback_render_atoms(self, context):
        blpanel = bpy.context.scene.blpanel
        print('Callback_render_atoms')
        render_atoms(blpanel.collection_list, blpanel.output_image)
    
    collection_list: EnumProperty(
        name="Collection",
        description="Collection",
        items=Callback_collection_list)
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
    output: StringProperty(
        name = "Output", default='batoms-output.xyz',
        description = "Output file")
    atoms_str: StringProperty(
        name = "Formula", default='H2O',
        description = "atoms_str")
    atoms_name: StringProperty(
        name = "Name", default='h2o',
        description = "name")
    output_image: StringProperty(
        name = "Output image", default='batoms.png',
        description = "output render image", update = Callback_render_atoms)
    

class AddMolecule(Operator):
    bl_idname = "batoms.add_molecule"
    bl_label = "Add molecule"
    bl_description = ("Add molecule")
    def execute(self, context):
        #print(molecule_str)
        blpanel = context.scene.blpanel
        atoms = molecule(blpanel.atoms_str)
        Batoms(label = blpanel.atoms_name, atoms = atoms, model_type = blpanel.model_type_add)
        return {'FINISHED'}
class AddBulk(Operator):
    bl_idname = "batoms.add_bulk"
    bl_label = "Add bulk"
    bl_description = ("Add bulk")
    def execute(self, context):
        blpanel = context.scene.blpanel
        atoms = bulk(blpanel.atoms_str)
        import_batoms(atoms, name = blpanel.atoms_name, model_type = blpanel.model_type_add, search_pbc_atoms={}, show_unit_cell=True)
        return {'FINISHED'}     
class AddAtoms(Operator):
    bl_idname = "batoms.add_atoms"
    bl_label = "Add atoms"
    bl_description = ("Add atoms")
    def execute(self, context):
        blpanel = context.scene.blpanel
        atoms = Atoms(blpanel.atoms_str)
        import_batoms(atoms, name = blpanel.atoms_name, model_type = blpanel.model_type_add)
        return {'FINISHED'}     

def modify_model_type(collection_name, model_type):
    batoms = Batoms(label = collection_name)
    batoms.model_type = model_type
def modify_scale(collection_name, scale):
    batoms = Batoms(label = collection_name)
    for species, ba in batoms.batoms.items():
        ba.scale = scale
def render_atoms(collection_name, output_image = 'bout.png'):
    batoms = Batoms(label = collection_name)
    batoms.render.run(output = output_image)