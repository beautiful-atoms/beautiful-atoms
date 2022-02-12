from ase.build import molecule, bulk
from ase import Atoms
import bpy
from bpy.types import (Panel,
                       Operator,
                       )
from bpy.props import (StringProperty,
                       IntProperty,
                       IntVectorProperty,
                       FloatProperty,
                       )
from batoms.butils import get_selected_batoms
from batoms import Batoms

# The panel.
class ASE_PT_prepare(Panel):
    bl_label       = "ASE"
    bl_space_type  = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "ASE_PT_Tools"

    def draw(self, context):
        layout = self.layout
        asepanel = context.scene.asepanel

        layout.prop(asepanel, "formula")
        layout.prop(asepanel, "label")
        layout.operator("batoms.add_molecule")
        layout.operator("batoms.add_bulk")
        layout.operator("batoms.add_atoms")
        #
        layout.label(text="Surface:")
        layout.prop(asepanel, "indices")
        layout.prop(asepanel, "nlayer")
        layout.prop(asepanel, "vacuum")
        layout.prop(asepanel, "termination")
        layout.operator("batoms.add_surface")



class ASEProperties(bpy.types.PropertyGroup):
    formula: StringProperty(
        name = "Formula", default='H2O',
        description = "formula")
    label: StringProperty(
        name = "Label", default='h2o',
        description = "Label")
    indices: IntVectorProperty(
        name="Miller indices", size = 3, default=(0, 0, 1),
        description = "Miller indices for the plane")
    vacuum: FloatProperty(
        name="vacuum", default=5.0,
        description = "vacuum")
    nlayer: IntProperty(name="layers", default=4)
    termination: StringProperty(
        name = "termination", default='',
        description = "termination")
    

class AddMolecule(Operator):
    bl_idname = "batoms.add_molecule"
    bl_label = "Add molecule"
    bl_description = ("Add molecule")
    def execute(self, context):
        asepanel = context.scene.asepanel
        atoms = molecule(asepanel.formula)
        Batoms(label = asepanel.label, from_ase = atoms)
        return {'FINISHED'}

class AddBulk(Operator):
    bl_idname = "batoms.add_bulk"
    bl_label = "Add bulk"
    bl_description = ("Add bulk")
    def execute(self, context):
        asepanel = context.scene.asepanel
        atoms = bulk(asepanel.formula)
        Batoms(label = asepanel.label, from_ase = atoms)
        return {'FINISHED'}  

class AddAtoms(Operator):
    bl_idname = "batoms.add_atoms"
    bl_label = "Add atoms"
    bl_description = ("Add atoms")
    def execute(self, context):
        asepanel = context.scene.asepanel
        atoms = Atoms(asepanel.formula)
        Batoms(label = asepanel.label, from_ase = atoms)
        return {'FINISHED'}
        
class AddSurface(Operator):
    @property
    def selected_batoms(self):
        return get_selected_batoms()
    bl_idname = "batoms.add_surface"
    bl_label = "Add surface"
    bl_description = ("Add surface")
    def execute(self, context):
        from ase.build import surface
        from ase.build.surfaces_with_termination import surfaces_with_termination
        asepanel = context.scene.asepanel
        selected_batoms = self.selected_batoms
        if len(selected_batoms) != 1:
            raise Exception('Please select one structure')
        batoms = Batoms(selected_batoms[0])
        bulk = batoms.atoms
        if len(asepanel.termination) == 0:
            atoms = surface(bulk, asepanel.indices, asepanel.nlayer, asepanel.vacuum)
        else:
            atoms = surfaces_with_termination(bulk, asepanel.indices, 
                                    asepanel.nlayer, asepanel.vacuum,
                                    termination=asepanel.termination)
        batoms.hide = True
        label = batoms.label + ''.join(str(i) for i in asepanel.indices)
        Batoms(label = label, from_ase = atoms, movie = True)
        return {'FINISHED'}