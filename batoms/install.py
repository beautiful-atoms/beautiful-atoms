"""

0) update batoms
1) set custom file folder.
2) set custom windows

"""

import bpy
from bpy.types import AddonPreferences
from bpy.props import (
    BoolProperty,
)
import pathlib
import subprocess
import sys
import os

account_name = "superstar54"
repo_name = "beautiful-atoms"
DEFAULT_PLUGIN_NAME = "batoms"

repo_git = f"https://github.com/{account_name}/{repo_name}.git"

def run(cmds):
    try:
        p = subprocess.Popen(cmds,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)
        stdout, stderr = p.communicate()
        if p.returncode == 0:
            return True
    except (OSError, Exception) as exception:
        return False
    return False

def has_git():
    return run(['git','--version'])

def gitclone(workdir=".", version="main", url=repo_git):
    """Make a git clone to the directory
    version can be a branch name or tag name
    """
    import shutil
    import batoms
    workdir = batoms.__path__[0]
    # workdir = pathlib.Path().resolve()
    addon_dir = os.path.dirname(workdir)
    clone_into = os.path.join(addon_dir, repo_name)
    print(clone_into)
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
    if run(commands):
        print(f"Cloned repo into directory {clone_into}")
    else:
        print("Faild to clone")
    #
    src = os.path.join(clone_into, "batoms")
    dst = os.path.dirname(workdir)
    if os.path.exists(workdir) and os.path.isdir(workdir):
        shutil.rmtree(workdir)
    print(src, dst)
    shutil.move(src, dst)

class BatomsUpdateButton(bpy.types.Operator):
    """Update Batoms"""
    bl_idname = "batoms.update"
    bl_label = "Update Batoms"
    bl_description = "Update to the latest development version."

    def execute(self, context):
        if not has_git():
            self.report({"ERROR"}, "Please install Git first.")
            print("Please install Git first.")
            return {"CANCELLED"}
        gitclone()
        return {"FINISHED"}

classes = [BatomsUpdateButton]

def register_class():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)

def unregister_class():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

if __name__ == "__main__":
    print(has_git())
    gitclone()
    