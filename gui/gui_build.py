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
from batoms.butils import get_selected_batoms, get_selected_vertices, \
                        get_selected_objects
from batoms import Batoms, Batom

# The panel.
class Build_PT_prepare(Panel):
    bl_label       = "Build"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "BUILD_PT_Tools"

    def draw(self, context):
        layout = self.layout
        bupanel = context.scene.bupanel

        box = layout.box()
        col = box.column(align=True)
        col.label(text="Add structure")
        col.prop(bupanel, "formula")
        col.prop(bupanel, "label")
        col.operator("batoms.add_molecule")
        col.operator("batoms.add_bulk")
        col.operator("batoms.add_atoms")


class BuildProperties(bpy.types.PropertyGroup):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    @property
    def selected_batom(self):
        return get_selected_objects('batom')
    @property
    def selected_vertices(self):
        return get_selected_vertices()
    
    
    formula: StringProperty(
        name = "Formula", default='H2O',
        description = "formula")
    label: StringProperty(
        name = "Label", default='h2o',
        description = "Label")

class AddMolecule(Operator):
    bl_idname = "batoms.add_molecule"
    bl_label = "Add molecule"
    bl_description = ("Add molecule")
    def execute(self, context):
        #print(molecule_str)
        bupanel = context.scene.bupanel
        atoms = molecule(bupanel.formula)
        Batoms(label = bupanel.label, atoms = atoms)
        return {'FINISHED'}
class AddBulk(Operator):
    bl_idname = "batoms.add_bulk"
    bl_label = "Add bulk"
    bl_description = ("Add bulk")
    def execute(self, context):
        bupanel = context.scene.bupanel
        atoms = bulk(bupanel.formula)
        Batoms(label = bupanel.label, atoms = atoms)
        return {'FINISHED'}     
class AddAtoms(Operator):
    bl_idname = "batoms.add_atoms"
    bl_label = "Add atoms"
    bl_description = ("Add atoms")
    def execute(self, context):
        bupanel = context.scene.bupanel
        atoms = Atoms(bupanel.formula)
        Batoms(label = bupanel.label, atoms = atoms)
        return {'FINISHED'}
