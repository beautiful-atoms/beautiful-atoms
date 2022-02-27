#!/usr/bin/env python3
import os
from os.path import expanduser, expandvars
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from xml.dom import NotFoundErr


minimal_version = 3.0
factory_python_directory = "_old_python"
repo_name = "beautiful-atoms"
repo_git = "https://github.com/superstar54/beautiful-atoms.git"
install_script_path = (
    "https://raw.githubusercontent.com/superstar54/beautiful-atoms/main/install.py"
)
plugin_name = "batoms"
BLENDERPY_ENABLE_PLUGIN = f"""
import bpy
import addon_utils
addon_utils.enable('{plugin_name}', default_set=True)
bpy.ops.wm.save_userpref()
print('Successfully enabled plugin {plugin_name}')
"""

BLENDERPY_DISABLE_PLUGIN = f"""
import bpy
import addon_utils
addon_utils.disable('{plugin_name}', default_set=True, handle_error=None)
bpy.ops.wm.save_userpref()
print('Disabled plugin {plugin_name}')
"""

BLENDERPY_TEST_PLUGIN = f"""
from batoms import Batoms
b = Batoms('O', ['O'], [[0, 0, 0]])
print('Test plugin import successful.')
"""

BLENDERPY_TEST_UNINSTALL = f"""
try:
    from batoms import Batoms
    raise Exception('batoms plugin still exists.')
except ImportError:
    print('batoms cleanly uninstalled.')
"""


def _get_default_locations(os_name, version=minimal_version):
    os_name = os_name.lower()
    if os_name not in ["windows", "macos", "linux"]:
        raise ValueError(f"{os_name} is not valid.")
    default_locations = {
        "windows": ["%PROGRAMFILES%/Blender Foundation/Blender {version}/{version}"],
        "macos": [
            "/Applications/Blender.app/Contents/Resources/{version}",
            "~/Applications/Blender.app/Contents/Resources/{version}",
        ],
    }
    if os_name == "linux":
        raise NotImplementedError(
            (
                "Version search in linux is not supported."
                " Please provide the full path to the blender installation location."
            )
        )
    matches = []
    for l in default_locations[os_name]:
        true_location = Path(
            expandvars(expanduser(l.format(version=version)))
        ).absolute()
        if true_location.is_dir():
            matches.append(true_location)
    if len(matches) > 1:
        print(f"Multiple blender installations of version {version} exists:")
        for i, m in enumerate(matches):
            print(f"{i}: {m.as_posix()}")
        choice = int(input("Choose one (default 0): ") or "0")
        match = matches[choice]
    elif len(matches) == 1:
        match = matches[0]
    else:
        match = None
    if match is None:
        raise FileNotFoundError(
            (
                f"No corresponding blender of version {version} on {os_name} found. "
                "Please specify the full path to the blender installation location."
            )
        )
    return match


def _get_blender_bin(os_name, blender_root):
    """Find the system-dependent blender binary file.
    blender_root is the path looks like <root>/3.0/
    """
    blender_root = Path(blender_root)
    if os_name == "windows":
        blender_bin = blender_root.parent / "blender.exe"
    elif os_name == "linux":
        blender_bin = blender_root.parent / "blender"
    elif os_name == "macos":
        # Blender.app/Contents/MacOS/Blender
        # Blender.app/Contents/Resources/3.0
        blender_bin = blender_root.parents[1] / "MacOS" / "Blender"
    else:
        raise NotImplementedError(f"Blender not supported for system {os_name}")

    if not blender_bin.exists():
        raise FileNotFoundError(
            f"Cannot find blender binary at {blender_bin.as_posix()}"
        )
    return blender_bin


def _get_os_name():
    p_ = sys.platform
    if p_ in ["win32"]:
        return "windows"
    elif p_ in ["linux", "linux2"]:
        return "linux"
    elif p_ in ["darwin"]:
        return "macos"
    else:
        raise NotImplementedError(f"Unsupported platform {p_}")


def is_conda():
    return "CONDA_PREFIX" in os.environ.keys()


def _get_conda_variables():
    results = {
        key: os.environ.get(key, "")
        for key in ("CONDA_PREFIX", "CONDA_PYTHON_EXE", "CONDA_DEFAULT_ENV")
    }
    return results


def _is_empty_dir(p):
    """Determin if path is empty"""
    p = Path(p)
    try:
        next(p.rglob("*"))
        stat = False
    except StopIteration:
        stat = True
    return stat


def _blender_enable_plugin(blender_bin):
    """Use blender's internal libary to enable plugin (and save as user script)"""
    commands = [blender_bin, "-b", "--python-expr", BLENDERPY_ENABLE_PLUGIN]
    proc = subprocess.run(commands)
    if proc.returncode != 0:
        raise RuntimeError("Enable plugin failed with error {proc.stderr}")
    return


def _blender_disable_plugin(blender_bin):
    """Use blender's internal libary to disable plugin (and save as user script)"""
    commands = [blender_bin, "-b", "--python-expr", BLENDERPY_DISABLE_PLUGIN]
    proc = subprocess.run(commands)
    if proc.returncode != 0:
        raise RuntimeError("Disable plugin failed with error {proc.stderr}")
    return


def gitclone(workdir=".", version="main", url=repo_git):
    """Make a git clone to the directory"""
    workdir = Path(expanduser(expandvars(workdir)))
    clone_into = workdir / repo_name
    os.makedirs(clone_into, exist_ok=True)
    commands = [
        "git",
        "clone",
        "--depth",
        "1",
        "-b",
        f"{version}",
        f"{url}",
        clone_into.as_posix(),
    ]
    proc = subprocess.run(commands)
    if proc.returncode == 0:
        print(f"Cloned repo into directory {clone_into.as_posix()}")
        return clone_into
    else:
        raise RuntimeError(f"Clone repo failed with error {proc.stderr}")


def install(
    blender_root,
    blender_bin,
    repo_path,
    factory_python_target=factory_python_directory,
    plugin_path_target=f"scripts/addons_contrib/{plugin_name}",
):
    """Copy the contents inside plugin_path to the target_directory"""
    blender_root = Path(blender_root)
    conda_env_file = repo_path / "env.yml"
    plugin_path_source = repo_path / plugin_name
    plugin_path_target = blender_root / plugin_path_target
    factory_python_source = blender_root / "python"
    factory_python_target = blender_root / factory_python_target

    # Check conda environment status
    conda_vars = _get_conda_variables()
    print(conda_vars)
    # Give a warning about conda env
    if conda_vars["CONDA_DEFAULT_ENV"] in ["base"]:
        print(
            (
                "Seems you're installing into the base environment. "
                "Installing batoms dependencies may interrupt your base environment. "
            )
        )
        choice = str(input("Continue? [y/N]") or "N").lower().startswith("y")
        if not choice:
            # TODO: update installation instruction
            print(
                "Abort. Please check the installation manual about how to activate an additional conda environment."
            )
            sys.exit(1)

    # Install from the env.yaml
    print("Updating conda environment")
    commands = [
        "conda",
        "env",
        "update",
        "-n",
        conda_vars["CONDA_DEFAULT_ENV"],
        "--file",
        conda_env_file.as_posix(),
    ]
    proc = subprocess.run(commands)
    if proc.returncode != 0:
        raise RuntimeError(f"Error updating conda env. Error is {proc.stderr}")

    # Symlink the conda env python --> `blender_root` / "python"
    # scenarios:
    # target _python is not empty:
    #       if _python is symlink --> unlink _python
    #       if _python is real path --> choose if overwrite
    # source python:
    #       if python is symlink --> link origin to _python if it does not exist (not likely from factory)
    #       if python is real path --> move to _python
    # If the "_python" path exists and not empty, only replace it if is symlink

    # Target part
    overwrite = False
    if factory_python_target.is_dir():
        if _is_empty_dir(factory_python_target):
            os.rmdir(factory_python_target)
        else:
            if factory_python_target.is_symlink():
                os.unlink(factory_python_target)
            else:
                print(
                    (
                        f"Target path {factory_python_target.as_posix()} is not empty. "
                        "This means you may have used the installation script already. "
                    )
                )
                try:
                    old_py = next(factory_python_target.glob("bin/python*"))
                except StopIteration:
                    print("Old python binary not found, may be broken.")
                    old_py = None
                    overwrite = True
                if old_py:
                    proc = subprocess.run([old_py, "-V"])
                    if proc.returncode == 0:
                        print("Old python binary is working. Will not overwrite.")
                        overwrite = False
                    else:
                        overwrite = True
    else:
        if factory_python_target.is_file():
            os.unlink(factory_python_target)

    # Source part
    if factory_python_target.exists():
        if overwrite:
            # Should not happen
            raise RuntimeError(
                f"Cannot overwrite {factory_python_target.as_posix()}. Check permission?"
            )
        else:
            pass
    else:
        if factory_python_source.is_symlink():
            origin = factory_python_source.readlink()
            os.unlink(factory_python_source)
            os.symlink(origin, factory_python_target)
        else:
            shutil.move(factory_python_source, factory_python_target)
        # finally, link the conda prefix of current environment
        conda_prefix = Path(conda_vars["CONDA_PREFIX"]).resolve()
        os.symlink(conda_prefix, factory_python_source)

    # Shall we overwrite the target path?
    if not _is_empty_dir(plugin_path_target):
        print(
            f"Target plugin installtion directory {plugin_path_target.as_posix()} is not empty."
        )
        choice = str(input("Overwrite? [y/N]") or "N").lower().startswith("y")
        if not choice:
            print("Abort.")
            sys.exit(1)
        else:
            if plugin_path_target.is_symlink():
                os.unlink(plugin_path_target)
            else:
                shutil.rmtree(plugin_path_target)
    shutil.copytree(plugin_path_source, plugin_path_target)
    print(f"Plugin copied to {plugin_path_target.as_posix()}.")

    _blender_enable_plugin(blender_bin)
    return


def uninstall(
    blender_root,
    blender_bin,
    factory_python_target=factory_python_directory,
    plugin_path_target=f"scripts/addons_contrib/{plugin_name}",
):
    """Remove the plugin from target_directory and restore"""
    blender_root = Path(blender_root)
    plugin_path_target = blender_root / plugin_path_target
    factory_python_source = blender_root / "python"
    factory_python_target = blender_root / factory_python_target
    _blender_disable_plugin(blender_bin)
    if plugin_path_target.is_symlink():
        os.unlink(plugin_path_target)
    elif plugin_path_target.is_dir():
        shutil.rmtree(plugin_path_target)
    else:
        print(
            f"Plugin directory {plugin_path_target.as_posix()} does not exist. Ignore."
        )

    # _old_python not found, ignore
    if not factory_python_target.exists():
        print(
            f"Backup of factory blender python path {factory_python_target.as_posix()} does not exist. Ignore"
        )
        return
    else:
        if factory_python_source.is_dir():
            if not _is_empty_dir(factory_python_source):
                if factory_python_source.is_symlink():
                    os.unlink(factory_python_source)
                elif factory_python_source.is_dir():
                    print(f"Current blender python path is not a symlink.")
                    overwrite = (
                        str(input("Overwrite? [y/N]") or "N").lower().startswith("y")
                    )
                    if overwrite:
                        shutil.rmtree(factory_python_source)
                    else:
                        print("Ignored.")
                        return
                else:
                    pass
            else:
                os.rmdir(factory_python_source)

        if factory_python_target.is_symlink():
            origin = factory_python_target.readlink()
            os.symlink(origin, factory_python_source)
        else:
            shutil.move(factory_python_target, factory_python_source)
    return


def test_plugin(blender_bin):
    commands = [
        blender_bin,
        "-b",
        "--python-exit-code",
        "1",
        "--python-expr",
        BLENDERPY_TEST_PLUGIN,
    ]
    proc = subprocess.run(commands)
    if proc.returncode != 0:
        raise RuntimeError("Plugin test failed with error {proc.stderr}")
    return


def test_uninstall(blender_bin):
    commands = [
        blender_bin,
        "-b",
        "--python-exit-code",
        "1",
        "--python-expr",
        BLENDERPY_TEST_UNINSTALL,
    ]
    proc = subprocess.run(commands)
    if proc.returncode != 0:
        raise RuntimeError("Uninstall test failed with error {proc.stderr}")
    return


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "blender_root",
        nargs="?",
        help=(
            "Path to the root of blender installation. "
            "If not provided, infer from os-dependent directories."
        ),
    )
    parser.add_argument(
        "-v", "--version", default=minimal_version, help="Blender major version number"
    )
    parser.add_argument(
        "-b",
        "--plugin-version",
        default="main",
        help="Plugin version or git hash tag. Default to fetch main",
    )
    parser.add_argument(
        "-p",
        "--local-repo-path",
        help="Path of local path for the repo, such as use git cloned.",
    )
    parser.add_argument(
        "--uninstall", action="store_true", help="Uninstall plugin in blender_root"
    )
    args = parser.parse_args()
    print(args)
    os_name = _get_os_name()
    # TODO: version compare!
    # Use blender_root supplied to commandline if exists, otherwise try to guess
    if args.blender_root:
        true_blender_root = Path(expanduser(expandvars(args.blender_root)))
    else:
        true_blender_root = _get_default_locations(os_name, version=args.version)
    true_blender_bin = _get_blender_bin(os_name, true_blender_root)
    print(f"Found blender binary at {true_blender_bin.as_posix()}")
    print(f"Choose blender directory at {true_blender_root.as_posix()}")
    if args.uninstall:
        uninstall(true_blender_root, true_blender_bin)
        test_uninstall(true_blender_bin)
    else:
        if not is_conda():
            print(
                "The installation script should be run inside a conda environment. Abort."
            )
            sys.exit(1)

        with tempfile.TemporaryDirectory() as workdir:
            if hasattr(args, "local_repo_path"):
                repo_path = Path(expanduser(expandvars(args.local_repo_path)))
            else:
                repo_path = gitclone(workdir, version=args.plugin_version, url=repo_git)
            install(true_blender_root, true_blender_bin, repo_path)
            test_plugin(true_blender_bin)


if __name__ == "__main__":
    main()
