import bpy
from ase import Atom, Atoms



def object_mode():
    for object in bpy.data.objects:
            if object.mode == 'EDIT':
                bpy.ops.object.mode_set(mode = 'OBJECT')

def read_batoms_collection_list():
    """
    Read all batoms collection 
    """
    items = [col.name for col in bpy.data.collections if col.batoms.is_batoms]
    return items
def read_atoms_list(coll):
    '''   
    '''
    coll_atom_kinds = [coll for coll in coll.children if 'atoms' in coll.name][0]
    elements = []
    for obj in coll_atom_kinds.all_objects:
        ele = obj.name.split('_')[2]
        elements.append(ele)
    return elements
        

def read_atoms_select():
    '''   
    '''
    from ase import Atoms, Atom
    atoms = Atoms()
    cell_vertexs = []
    for obj in bpy.context.selected_objects:
        if "BOND" in obj.name.upper():
            continue
        if obj.type not in {'MESH', 'SURFACE', 'META'}:
            continue
        name = ""
        if 'atom_kind_' == obj.name[0:10]:
            print(obj.name)
            ind = obj.name.index('atom_kind_')
            ele = obj.name[ind + 10:].split('_')[0]
            if len(obj.children) != 0:
                for vertex in obj.data.vertices:
                    location = obj.matrix_world @ vertex.co
                    atoms.append(Atom(ele, location))
            else:
                if not obj.parent:
                    location = obj.location
                    atoms.append(Atom(ele, location))
        # cell
        if 'point_cell' == obj.name[0:10]:
            print(obj.name)
            if len(obj.children) != 0:
                for vertex in obj.data.vertices:
                    location = obj.matrix_world @ vertex.co
                    # print(location)
                    cell_vertexs.append(location)
    # print(atoms)
    # print(cell_vertexs)
    if cell_vertexs:
        cell = [cell_vertexs[4], cell_vertexs[2], cell_vertexs[1]]
        atoms.cell = cell
        atoms.pbc = [True, True, True]
    # self.atoms = atoms
    return atoms

def remove_collection(name):
        collection = bpy.data.collections.get(name)
        for obj in collection.all_objects:
            bpy.data.objects.remove(obj, do_unlink=True)
        for coll in collection.children:
            bpy.data.collections.remove(coll)
        bpy.data.collections.remove(collection)
        
def clean_objects():
    for item in bpy.data.objects:
        bpy.data.objects.remove(item)

def removeAll():
    #types =  ['MESH', 'CURVE', 'SURFACE', 'META', 'FONT', 'ARMATURE', 'LATTICE', 'EMPTY', 'GPENCIL', 'CAMERA', 'LIGHT', 'SPEAKER', 'LIGHT_PROBE']
    for mesh in bpy.data.meshes:
        bpy.data.meshes.remove(mesh)
    for obj in bpy.data.objects:
        bpy.data.objects.remove(obj)
    for cam in bpy.data.cameras:
        bpy.data.cameras.remove(cam)
    for light in bpy.data.lights:
        bpy.data.lights.remove(light)
    for mat in bpy.data.materials:
        bpy.data.materials.remove(mat)
    for coll in bpy.data.collections:
        if coll.name == 'Collection': continue
        bpy.data.collections.remove(coll)

def lock_camera_to_view(switch):
    for area in bpy.context.screen.areas:
        if area.type == 'VIEW_3D':
            for space in area.spaces:
                if space.type == 'VIEW_3D':
                    space.lock_camera = switch