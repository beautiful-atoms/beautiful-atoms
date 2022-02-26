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

minimal_version = 3.0
factory_python_directory = "_old_python"
repo_name = "beautiful-atoms"
repo_git = "https://github.com/superstar54/beautiful-atoms.git"
install_script_path = (
    "https://raw.githubusercontent.com/superstar54/beautiful-atoms/main/install.py"
)


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
    else:
        match = matches[0]
    return match


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
    results = dict(
        os.environ.get(key, "")
        for key in ("CONDA_PREFIX", "CONDA_PYTHON_EXE", "CONDA_DEFAULT_ENV")
    )
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
    repo_path,
    factory_python_target="factory_python_directory",
    plugin_path_target="scripts/addons_contrib/batoms",
):
    """Copy the contents inside plugin_path to the target_directory"""
    blender_root = Path(blender_root)
    conda_env_file = blender_root / "env.yml"
    plugin_path_source = repo_path / "batoms"
    plugin_path_target = blender_root / plugin_path_target
    # Shall we overwrite the target path?
    if not _is_empty_dir(plugin_path_target):
        print(f"Target plugin installtion directory {plugin_path_target.as_posix()} is not empty.")
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
    print(f"Plugin copied to  {plugin_path_target.as_posix()}.")

    pass


def uninstall(
    blender_root,
    plugin_name="batoms",
    original_python_directory="factory_python_directory",
    target_directory="scripts/addons_contrib/batoms",
):
    """Remove the plugin from target_directory"""
    pass


def test_plugin():
    pass


def test_uninstall():
    pass


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
        "--uninstall", action="store_true", help="Uninstall plugin in blender_root"
    )
    args = parser.parse_args()
    print(args)
    os_name = _get_os_name()
    true_blender_root = _get_default_locations(os_name, version=args.version)
    print(f"Choose blender directory at {true_blender_root.as_posix()}")
    if args.uninstall:
        uninstall(true_blender_root)
        test_uninstall()
    else:
        if not is_conda():
            raise RuntimeError(
                "The installation script should be run inside a conda environment."
            )

        with tempfile.TemporaryDirectory() as workdir:
            cloned_repo = gitclone(workdir, version=args.plugin_version, url=repo_git)
            install(true_blender_root, cloned_repo)
            test_plugin()


if __name__ == "__main__":
    main()
