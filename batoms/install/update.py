"""
update batoms
"""

import bpy
import subprocess
import os
import logging
# logger = logging.getLogger('batoms')
logger = logging.getLogger(__name__)

account_name = "beautiful-atoms"
repo_name = "beautiful-atoms"
DEFAULT_PLUGIN_NAME = "batoms"

repo_git = f"https://github.com/{account_name}/{repo_name}.git"

def subprocess_run(cmds):
    try:
        p = subprocess.Popen(cmds,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)
        stdout, stderr = p.communicate()
        if p.returncode == 0:
            return True
        else:
            logger.critical(stdout, stderr)
    except (OSError, Exception) as exception:
        logger.critical(exception)
        return False
    return False

def has_git():
    return subprocess_run(['git','--version'])

def gitclone(workdir=".", version="main", url=repo_git):
    """Make a git clone to the directory
    version can be a branch name or tag name
    """
    import shutil
    import pathlib
    batoms_dir = os.path.dirname(pathlib.Path(__file__).parent.resolve())
    addon_dir = os.path.dirname(batoms_dir)
    clone_into = os.path.join(addon_dir, repo_name)
    if os.path.exists(clone_into) and os.path.isdir(clone_into):
        shutil.rmtree(clone_into)
    commands = [
        "git",
        "clone",
        "--depth",
        "1",
        "-b",
        f"{version}",
        f"{url}",
        clone_into,
    ]
    flag = subprocess_run(commands)
    flag = True
    if flag:
        logger.info(f"Cloned repo into directory {clone_into}")
    else:
        logger.warning("Faild to clone")
    #
    src = os.path.join(clone_into, "batoms")
    if os.path.exists(batoms_dir) and os.path.isdir(batoms_dir):
        logger.info("Remove old Batoms folder {}".format(batoms_dir))
        shutil.rmtree(batoms_dir)
    logger.info("Move {} to {}".format(src, batoms_dir))
    shutil.move(src, batoms_dir)

class BatomsUpdateButton(bpy.types.Operator):
    """Update Batoms"""
    bl_idname = "batoms.update"
    bl_label = "Update Batoms"
    bl_description = "Update to the latest development version."

    def execute(self, context):
        if not has_git():
            self.report({"ERROR"}, "Please install Git first.")
            logger.warning("Please install Git first.")
            return {"CANCELLED"}
        gitclone()
        self.report({"INFO"}, "Update to the latest version successfully!")
        return {"FINISHED"}
