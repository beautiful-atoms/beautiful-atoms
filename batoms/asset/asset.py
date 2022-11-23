from pathlib import Path
import bpy
import pathlib
import os
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)



batoms_asset_dir = os.path.join(pathlib.Path(__file__).parent.resolve(), 'libraries')


def update_asset_preview():
    for coll in bpy.data.collections:
        if coll.asset_data is not None:
            coll.asset_generate_preview()

def save_asset(filepath, assets):
    bpy.data.libraries.write(filepath, assets)

def push_asset(library, file, assets):
    from pathlib import Path
    import shutil
    # save name of assets
    asset_names = [asset.name for asset in assets]
    # save to local path
    localpath = './tmp.blend'
    save_asset(localpath, assets)
    #
    path = bpy.context.preferences.filepaths.asset_libraries[library].path
    filepath = os.path.join(path, file)
    if not os.path.exists(filepath):
        emptypath = os.path.join(path, 'empty.blend')
        shutil.copyfile(emptypath, filepath)
    bpy.ops.wm.open_mainfile(filepath=filepath)
    localpath = Path(localpath)
    for name in asset_names:
        bpy.ops.wm.append(
            filepath=str(localpath / 'Collection' / name),
            directory=str(localpath / 'Collection'),
            filename=name,
            )
    bpy.ops.wm.save_mainfile(filepath=filepath, compress=True)
    bpy.ops.wm.open_mainfile(filepath='./tmp.blend')

def find_baotms_catalog_uuid(library, catalog_name, directory = None):
    """Find catalog by name.

    Args:
        catalog_name (str): Name of the catalog
        directory (str, optional): _description_. Defaults to None.

    Returns:
        str: uuid of the catalog.
    """
    if directory is None:
        directory = bpy.context.preferences.filepaths.asset_libraries[library].path
    logger.debug(directory)
    asset_cats = os.path.join(directory, "blender_assets.cats.txt")
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
    import time
    batoms.model_style = model_style
    batoms.coll.asset_mark()
    batoms.coll.asset_generate_preview()
    batoms.coll.asset_data.author = metadata["author"]
    batoms.coll.asset_data.description = metadata["description"]
    batoms.coll.asset_data.catalog_id = metadata["catalog_id"]
    for tag in metadata["tags"]:
        batoms.coll.asset_data.tags.new(name=tag)

    return batoms.coll
