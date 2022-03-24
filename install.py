#!/usr/bin/env python3
"""
Universal installation script for batoms using anaconda or miniconda environment
Installation can be as easy as following steps
```
conda create -n batoms git curl python
conda activate batoms
curl -O https://raw.githubusercontent.com/superstar54/beautiful-atoms/main/install.py
python install.py
```

Alternatively, check what `install.py` can do by using
```
python install.py --help
```
"""
import os
from os.path import expanduser, expandvars
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from distutils.version import LooseVersion

# TODO: allow version control
# TODO: windows privilege issue
# TODO: complete install tutorial about the env variables
# TODO: unit test workflow


DEFAULT_GITHUB_ACCOUNT = "superstar54"
DEFAULT_REPO_NAME = "beautiful-atoms"
DEFAULT_PLUGIN_NAME = "batoms"
DEFAULT_PLUGIN_PATH = f"scripts/addons_contrib/{DEFAULT_PLUGIN_NAME}"

MIN_BLENDER_VER = "3.0"
PY_PATH = "python"
PY_BACKUP_PATH = "_old_python"

BLENDERPY_ENABLE_PLUGIN = f"""
import bpy
import addon_utils
addon_utils.enable('{DEFAULT_PLUGIN_NAME}', default_set=True)
bpy.ops.wm.save_userpref()
print('Successfully enabled plugin {DEFAULT_PLUGIN_NAME}')
"""

BLENDERPY_DISABLE_PLUGIN = f"""
import bpy
import addon_utils
addon_utils.disable('{DEFAULT_PLUGIN_NAME}', default_set=True, handle_error=None)
bpy.ops.wm.save_userpref()
print('Disabled plugin {DEFAULT_PLUGIN_NAME}')
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

# The directory to move factory python from conda

repo_name = os.environ.get("GITHUB_REPO", DEFAULT_REPO_NAME)
account_name = os.environ.get("GITHUB_ACCOUNT", DEFAULT_GITHUB_ACCOUNT)

repo_git = f"https://github.com/{account_name}/{repo_name}.git"
install_script_path = (
    f"https://raw.githubusercontent.com/{account_name}/{repo_name}/main/install.py"
)


def _get_default_locations(os_name, version=MIN_BLENDER_VER):
    """Get system specific default install locations of blender"""
    os_name = os_name.lower()
    # Compare version
    if LooseVersion(str(version)) < LooseVersion(str(MIN_BLENDER_VER)):
        raise ValueError(
            f"Blender version {version} is not supported. Minimal requirement is {MINIMAL_BLENDER_VER}"
        )
    if os_name not in ["windows", "macos", "linux"]:
        raise ValueError(f"{os_name} is not valid.")
    default_locations = {
        "windows": ["%PROGRAMFILES%/Blender Foundation/Blender {version}/{version}"],
        "macos": [
            "/Applications/Blender.app/Contents/Resources/{version}",
            "$HOME/Applications/Blender.app/Contents/Resources/{version}",
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
    # Multiple blender installation may occur on macos
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


def _get_blender_bin(os_name, blender_bundle_root):
    """Find the system-dependent blender binary file.
    blender_bundle_root contains the distribution bundle
    and has the pattern <blender_root>/<version>/
    """
    blender_bundle_root = Path(blender_bundle_root)
    if os_name == "windows":
        blender_bin = blender_bundle_root.parent / "blender.exe"
    elif os_name == "linux":
        blender_bin = blender_bundle_root.parent / "blender"
    elif os_name == "macos":
        # Blender.app/Contents/MacOS/Blender
        # Blender.app/Contents/Resources/3.0
        blender_bin = blender_bundle_root.parents[1] / "MacOS" / "Blender"
    else:
        raise NotImplementedError(f"Blender not supported for system {os_name}")

    if not blender_bin.exists():
        raise FileNotFoundError(
            f"Cannot find blender binary at {blender_bin.as_posix()}"
        )
    return blender_bin


def _get_os_name():
    """Convient os name function"""
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
        for key in (
            "CONDA_PREFIX",
            "CONDA_PYTHON_EXE",
            "CONDA_DEFAULT_ENV",
            "CONDA_EXE",
        )
    }
    return results


def _is_empty_dir(p):
    """Determin if path has no child dirs"""
    p = Path(p)
    try:
        next(p.rglob("*"))
        stat = False
    except StopIteration:
        stat = True
    return stat


def _run_process(commands, shell=False, print_cmd=True, cwd="."):
    full_cmd = " ".join(commands)
    if print_cmd:
        print(" ".join(commands))
    proc = subprocess.run(commands, shell=shell, cwd=cwd)
    if proc.returncode == 0:
        return proc.returncode
    else:
        raise RuntimeError(f"Running {full_cmd} returned error code {proc.stderr}")


def _blender_enable_plugin(blender_bin):
    """Use blender's internal libary to enable plugin (and save as user script)"""
    blender_bin = str(blender_bin)
    commands = [blender_bin, "-b", "--python-expr", BLENDERPY_ENABLE_PLUGIN]
    _run_process(commands)
    return


def _blender_disable_plugin(blender_bin):
    """Use blender's internal libary to disable plugin (and save as user script)"""
    blender_bin = str(blender_bin)
    commands = [blender_bin, "-b", "--python-expr", BLENDERPY_DISABLE_PLUGIN]
    _run_process(commands)
    return


def _gitclone(workdir=".", version="main", url=repo_git):
    """Make a git clone to the directory
    version can be a branch name or tag name
    """
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
    _run_process(commands)
    print(f"Cloned repo into directory {clone_into.as_posix()}")

def _rename_dir(src, dst):
    pass


def _conda_update(conda_env_file, conda_vars, env_name=None):
    """Update conda environment using env file.
    If env_name is None, use default env
    if reinstall_numpy, use pip to reinstall numpy (windows only)
    """
    conda_env_file = Path(conda_env_file)
    if env_name is None:
        env_name = conda_vars["CONDA_DEFAULT_ENV"]

    if env_name in ["base"]:
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
        conda_vars["CONDA_EXE"],
        "env",
        "update",
        "-n",
        env_name,
        "--file",
        conda_env_file.as_posix(),
    ]
    _run_process(commands)
    print("Finished install conda packages.")

    # Extra steps (windows only), replace the numpy version with pip's numpy
    # since python on windows comes with wheel this should be relatively straighforward
    if _get_os_name() in [
        "windows",
    ]:
        print(
            (
                "You're running on Windows. "
                "We will try to use wheel provided by pip for Numpy "
                "to resolve DLL not found issue."
            )
        )
        commands = [
            "python",
            "-m",
            "pip",
            "uninstall",
            "-y",
            "numpy",
        ]
        _run_process(commands)
        print("Uninstalled conda-distributed Numpy.")

        commands = [
            "python",
            "-m",
            "pip",
            "install",
            "-y",
            "numpy",
        ]
        _run_process(commands)
        print("Reinstalled Numpy from pip wheel")
    return


def install(blender_root, blender_bin, repo_path):
    """Link current conda environment to blender's python root
    Copy batoms plugin under repo_path to plugin directory
    """
    blender_root = Path(blender_root)
    conda_env_file = repo_path / "env.yml"
    # plugin_path_source: dir of plugin in github repo
    plugin_path_source = repo_path / DEFAULT_PLUGIN_NAME
    # plugin_path_target: dir of plugin under blender root
    plugin_path_target = blender_root / DEFAULT_PLUGIN_PATH
    # factory_python_source: (presumably) python shipped with blender
    factory_python_source = blender_root / PY_PATH
    # factory_python_target: dir to move factory python
    factory_python_target = blender_root / PY_BACKUP_PATH

    # Check conda environment status
    conda_vars = _get_conda_variables()
    
    # Installation logic follows the COMPAS project
    # If the factory_python_target exists, restore factory python first
    #    detect if the target is in a good state
    # Move the factory_python_source to factory_python_target
    # Symlink conda prefix to factory_python_source
    # TODO: wrap the directory handling into another function

    if factory_python_target.exists():
        # Empty of symlink factory python might be mistakenly created
        # just delete them
        exception_corrupt = False
        if _is_empty_dir(factory_python_target) or factory_python_target.is_symlink():
            exception_corrupt = True
        else:
            print(f"Target path {factory_python_target.as_posix()} exists.")
            try:
                old_py = next(factory_python_target.glob("bin/python*"))
            except StopIteration:
                print("Factory python binary not found.")
                old_py = None
                exception_corrupt = True
            if old_py:
                proc = subprocess.run([old_py, "-V"])
                if proc.returncode != 0:
                    print(
                        (
                            f"Found factory python at {old_py.as_posix()} "
                            "but it's not working"
                        )
                    )
                    exception_corrupt = True
        # TODO: improve error msg
        if exception_corrupt:
            raise RuntimeError(
                (
                    f"Backup of Blender's factory python at {factory_python_target.as_posix()} "
                    "is corrupt. Please replace it with a copy."
                )
            )
        else:
            if factory_python_source.exists():
                if factory_python_source.is_symlink():
                    os.unlink(factory_python_source)
                elif factory_python_source.is_dir():
                    raise OSError(
                        f"{factory_python_source.as_posix()} is not symlink. Please backup its content and retry."
                    )
                else:
                    os.unlink(factory_python_source)
            os.rename(factory_python_target, factory_python_source)
            print(
                f"Renamed {factory_python_target.as_posix()} to {factory_python_source.as_posix()}"
            )

    # Step 2: rename soruce to target
    if factory_python_target.exists():
        raise RuntimeError(
            f"Something wrong. {factory_python_target.as_posix()} still exists."
        )
    if (not factory_python_source.is_dir()) or (factory_python_source.is_symlink()):
        raise RuntimeError(
            f"Something wrong. {factory_python_source.as_posix()} should be a real directory."
        )
    os.rename(factory_python_source, factory_python_target)
    print(
        f"Renamed {factory_python_source.as_posix()} to {factory_python_target.as_posix()}"
    )

    # Step 3: finally, link the conda prefix of current environment
    conda_prefix = Path(conda_vars["CONDA_PREFIX"]).resolve()
    # Should not happen but just in case
    if factory_python_source.is_symlink():
        os.unlink(factory_python_source)
    os.symlink(conda_prefix, factory_python_source)
    print(f"Created symlink {conda_prefix} --> {factory_python_source.as_posix()}")

    # Give a warning about conda env
    # TODO: allow direct install into another environment
    _conda_update(conda_env_file, conda_vars)

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
):
    """Remove the plugin from target_directory and restore"""
    blender_root = Path(blender_root)
    plugin_path_target = blender_root / DEFAULT_PLUGIN_PATH
    factory_python_source = blender_root / PY_PATH
    factory_python_target = blender_root / PY_BACKUP_PATH
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
    blender_bin = str(blender_bin)
    commands = [
        blender_bin,
        "-b",
        "--python-exit-code",
        "1",
        "--python-expr",
        BLENDERPY_TEST_PLUGIN,
    ]
    _run_process(commands)
    return


def test_uninstall(blender_bin):
    blender_bin = str(blender_bin)
    commands = [
        blender_bin,
        "-b",
        "--python-exit-code",
        "1",
        "--python-expr",
        BLENDERPY_TEST_UNINSTALL,
    ]
    _run_process(commands)
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
        "-v", "--version", default=MIN_BLENDER_VER, help="Blender major version number"
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
                repo_path = _gitclone(
                    workdir, version=args.plugin_version, url=repo_git
                )
            install(true_blender_root, true_blender_bin, repo_path)
            test_plugin(true_blender_bin)


if __name__ == "__main__":
    main()
