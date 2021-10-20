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
        col.label(text="Replace by:")
        col.prop(bapanel, "species")

        box = layout.box()
        col = box.column(align=True)
        col.operator("batoms.measurement")
        col.prop(bapanel, "measurement")

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
    @property
    def selected_batom(self):
        return get_selected_objects('batom')
    @property
    def selected_vertices(self):
        return get_selected_vertices()
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
    def Callback_replace_atoms(self, context):
        bapanel = bpy.context.scene.bapanel
        replace_atoms(self.selected_vertices, bapanel.species)
    
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
    measurement: StringProperty(
        name="value", default='',
        description = "measurement in Angstrom, degree")
    species: StringProperty(
        name="", default='O_1',
        description = "replaced by species", update = Callback_replace_atoms)
    atoms_str: StringProperty(
        name = "Formula", default='H2O',
        description = "atoms_str")
    atoms_name: StringProperty(
        name = "Label", default='h2o',
        description = "Label")
    filetype: StringProperty(
        name = "File type", default='xyz',
        description = "save batoms to file", update = Callback_export_atoms)


def modify_batoms_attr(selected_batoms, key, value):
    """
    """
    batoms_list = []
    for name in selected_batoms:
        batoms = Batoms(label = name)
        setattr(batoms, key, value)
        batoms_list.append(batoms)
    # Restor the selected state
    for batoms in batoms_list:
        batoms.select = True
def export_atoms(selected_batoms, filetype = 'xyz'):
    batoms_list = []
    for name in selected_batoms:
        batoms = Batoms(label = name)
        batoms.write('%s.%s'%(batoms.label, filetype))
        batoms_list.append(batoms)
    for batoms in batoms_list:
        batoms.select = True

def replace_atoms(selected_vertices, species):
    batom_list = []
    for name, index in selected_vertices:
        obj = bpy.data.objects[name]
        batoms = Batoms(label = obj.batom.label)
        batoms.replace(obj.batom.species, species, index)
        batom_list.append(name)
    for name in batom_list:
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.context.view_layer.objects.active = obj
            bpy.ops.object.mode_set(mode='EDIT')
def measurement():
    """
    Todo: only distance works.

    We need histroy of selections.

    
    bmesh.select_history only works for one objects.
    """
    import numpy as np
    from ase.geometry.geometry import get_distances, get_angles, get_dihedrals
    selected_vertices = get_selected_vertices()
    cell = None
    pbc = None
    batom_list = []
    positions = np.array([]).reshape(-1, 3)
    for name, index in selected_vertices:
        obj = bpy.data.objects[name]
        batom = Batom(label = name)
        positions = np.append(positions, batom.positions[index], axis = 0)
        batom_list.append(name)
    if len(positions) == 2:
        results = get_distances([positions[0]], 
                    [positions[1]], 
                    cell=cell, pbc=pbc)[1]
    elif len(positions) == 3:
        v12 = positions[0] - positions[1]
        v32 = positions[2] - positions[1]
        results =  get_angles([v12], [v32], cell=cell, pbc=pbc)
    elif len(positions) == 4:
        v0 = positions[1] - positions[0]
        v1 = positions[2] - positions[1]
        v2 = positions[3] - positions[2]
        results =  get_dihedrals([v0], [v1], [v2], cell=cell, pbc=pbc)
    else:
        return 'Not supported'
    for name in batom_list:
        if name in bpy.data.objects:
            obj = bpy.data.objects.get(name)
            bpy.context.view_layer.objects.active = obj
            bpy.ops.object.mode_set(mode='EDIT')
    results.shape = (-1,)
    results = [str(round(float(i), 2)) for i in results]
    results = ' '.join(results)
    return results

class MeasureButton(Operator):
    bl_idname = "batoms.measurement"
    bl_label = "Measure"
    bl_description = "Measure distance, angle and dihedra angle"

    def execute(self, context):
        bapanel = context.scene.bapanel
        result = measurement()
        bapanel.measurement = result
        return {'FINISHED'}

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
        Batoms(label = bapanel.atoms_name, atoms = atoms)
        return {'FINISHED'}     
class AddAtoms(Operator):
    bl_idname = "batoms.add_atoms"
    bl_label = "Add atoms"
    bl_description = ("Add atoms")
    def execute(self, context):
        bapanel = context.scene.bapanel
        atoms = Atoms(bapanel.atoms_str)
        Batoms(label = bapanel.atoms_name, atoms = atoms)
        return {'FINISHED'}
