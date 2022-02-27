from cProfile import label
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
from batoms.utils.butils import get_selected_batoms
from batoms import Batoms

# The panel.


class ASE_PT_prepare(Panel):
    bl_label = "ASE"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_category = "Batoms"
    bl_idname = "ASE_PT_Tools"

    def draw(self, context):
        layout = self.layout
        asepanel = context.scene.asepanel

        layout.prop(asepanel, "formula")
        layout.prop(asepanel, "label")
        op = layout.operator("batoms.add_molecule")
        op.formula = asepanel.formula
        op.label = asepanel.label
        op = layout.operator("batoms.add_bulk")
        op.formula = asepanel.formula
        op.label = asepanel.label
        op = layout.operator("batoms.add_atoms")
        op.formula = asepanel.formula
        op.label = asepanel.label
        #
        layout.label(text="Surface:")
        layout.prop(asepanel, "indices")
        layout.prop(asepanel, "nlayer")
        layout.prop(asepanel, "vacuum")
        layout.prop(asepanel, "termination")
        op = layout.operator("batoms.add_surface")
        op.label = asepanel.label
        op.indices = asepanel.indices
        op.size = asepanel.size
        op.vacuum = asepanel.vacuum
        op.nlayer = asepanel.nlayer
        op.termination = asepanel.termination


class ASEProperties(bpy.types.PropertyGroup):
    formula: StringProperty(
        name="Formula", default='H2O',
        description="formula")
    label: StringProperty(
        name="Label", default='h2o',
        description="Label")
    indices: IntVectorProperty(
        name="Miller indices", size=3, default=(1, 1, 1),
        min = -20, max = 20,
        description="Miller indices for the surface")
    size: IntVectorProperty(
        name="Size", size=3, default=(1, 1, 1),
        min = 1, max = 20,
        description="System size in units of the minimal unit cell.")
    vacuum: FloatProperty(
        name="vacuum", default=5.0,
        min = 0, max = 100,
        description="vacuum")
    nlayer: IntProperty(name="layers", min = 1, max = 50, default=4)
    termination: StringProperty(
        name="termination", default='',
        description="termination")

