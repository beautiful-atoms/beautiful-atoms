from pathlib import Path
import bpy
import pathlib
import os

batoms_asset_dir = os.path.join(pathlib.Path(__file__).parent.resolve(), 'libraries')


def push_asset(library, file, formula = 'Pt'):
    path = bpy.context.preferences.filepaths.asset_libraries[library].path
    f = os.path.join(path, file)
    bpy.ops.wm.open_mainfile(filepath=f)
    bpy.ops.batoms.bulk_add(label = 'pt', formula = formula, cubic = True)
    bpy.ops.wm.save_mainfile()


def run_blender(script):
    blender_cmd = 'blender'
    # root = os.path.normpath(os.path.dirname(__file__))
    # script = os.path.join(root, 'script-api.py')
    cmd = blender_cmd + ' -b ' + ' -P ' + script
    errcode = os.system(cmd)
    if errcode != 0:
        raise OSError('Command ' + cmd +
                      ' failed with error code %d' % errcode)


def find_baotms_catalog_uuid(catalog_name, directory = None):
    """Find catalog by name.

    Args:
        catalog_name (str): Name of the catalog
        directory (str, optional): _description_. Defaults to None.

    Returns:
        str: uuid of the catalog.
    """
    if directory is None:
        batoms_asset_dir = directory
    asset_cats = os.path.join(batoms_asset_dir, "blender_assets.cats.txt")
    with open(asset_cats) as f:
        for line in f.readlines():
            if line.startswith(("#", "VERSION", "\n")):
                continue
            # uuid:catalog_tree:catalog_name'+'\n'
            name = line.rstrip().split(":")[2]
            if catalog_name == name:
                uuid = line.split(":")[0]
                return uuid

def create_asset(batoms, model_style = 0, metadata = {}):
    """Create a asset using batoms.

    Args:
        batoms (Batoms): _description_
        metadata (dict): _description_
    """
    batoms.model_style = model_style
    # hide instancer
    for obj in batoms.coll.all_objects[:]:
        if obj.batoms.type == "INSTANCER":
            obj.scale = [0.001, 0.001, 0.001]
    batoms.coll.asset_mark()
    batoms.coll.asset_generate_preview()
    batoms.coll.asset_data.author = metadata["author"]
    batoms.coll.asset_data.description = metadata["description"]
    batoms.coll.asset_data.catalog_id = metadata["catalog_id"]
    for tag in metadata["tags"]:
        batoms.coll.asset_data.tags.new(name=tag)
