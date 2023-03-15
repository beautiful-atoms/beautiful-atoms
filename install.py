#!/usr/bin/env python3
"""
Universal installation script for batoms using anaconda or miniconda environment
Installation can be as easy as following steps
```
conda create -n batoms git curl python
conda activate batoms
curl -O https://raw.githubusercontent.com/beautiful-atoms/beautiful-atoms/main/install.py
python install.py
```

Alternatively, check what `install.py` can do by using
```
python install.py --help
```
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from os.path import expanduser, expandvars

try:
    from packaging.version import Version
except ImportError:
    from distutils.version import LooseVersion as Version
# TODO: allow version control
# TODO: windows privilege issue
# TODO: complete install tutorial about the env variables
# TODO: unit test workflow
# TODO: python && numpy version match
# sys.version --> python version
# numpy.__version__ --> numpy

################################################################################
# Section 1: default variables
################################################################################

DEFAULT_GITHUB_ACCOUNT = "beautiful-atoms"
DEFAULT_REPO_NAME = "beautiful-atoms"
DEFAULT_PLUGIN_NAME = "batoms"
DEFAULT_PLUGIN_PATH = f"scripts/addons_contrib/{DEFAULT_PLUGIN_NAME}"
DEFAULT_BLENDER_PY_VER = "3.10.2"
DEFAULT_BLENDER_NUMPY_VER = "1.22.0"

ALLOWED_BLENDER_VERSIONS = ["3.4", "3.5"]
DEPRECATED_BLENDER_VERSIONS = ["3.0", "3.1", "3.2"]

PY_PATH = "python"
PY_BACKUP_PATH = "_old_python"
BLENDER_BACKUP_PATH = "blender.origin"  # Backup of the factory blender binary

REPO_NAME = os.environ.get("GITHUB_REPO", DEFAULT_REPO_NAME)
ACCOUNT_NAME = os.environ.get("GITHUB_ACCOUNT", DEFAULT_GITHUB_ACCOUNT)

REPO_GIT = f"https://github.com/{ACCOUNT_NAME}/{REPO_NAME}.git"
INSTALL_SCRIPT_PATH = (
    f"https://raw.githubusercontent.com/{ACCOUNT_NAME}/{REPO_NAME}/main/install.py"
)

################################################################################
# Section 2: temporary file script templates (python or yaml)
################################################################################

# Adding an embedded YAML file for better portability
ENV_YAML = """
name: beautiful_atoms
channels:
  - conda-forge
  - anaconda
  - defaults
dependencies:
  - git
  - curl
  - python={blender_py_ver}
  - numpy={blender_numpy_ver}
  - pymatgen>=2020.12
  - ase>=3.21.0
  - openbabel>=3.1.1
  - scikit-image
  - scipy>=1.6.0
  - spglib>=2.0.0
  - pip
  - pip:
    - batoms-api
"""

BLENDER_CHK_VERSION = """
import sys
import numpy
print("Python Version: ", sys.version)
print("Numpy Version: ", numpy.__version__)
"""

BLENDERPY_ENABLE_PLUGIN = f"""
print('start')
import bpy
print('bpy can import')
import addon_utils
print('add on can import')
print('Searching addons from: ', addon_utils.paths())
print(addon_utils.check('batoms'))
addon_utils.enable('{DEFAULT_PLUGIN_NAME}', default_set=True)
print(addon_utils.check('batoms'))
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

BLENDERPY_TEST_PLUGIN = """
from batoms import Batoms
b = Batoms('O', ['O'], [[0, 0, 0]])
print('Test plugin import successful.')
"""

BLENDERPY_TEST_UNINSTALL = """
try:
    from batoms import Batoms
    raise Exception('batoms plugin still exists.')
except ImportError:
    print('batoms cleanly uninstalled.')
"""

BLENDERPY_SETTING_PREFERENCES = """
print('start setting preferences')
import bpy
bpy.ops.batoms.use_batoms_preference()
print('Successfully setting and preference.')
"""

BLENDERPY_SETTING_STARTUP = """
print('start setting startup')
import bpy
bpy.ops.batoms.use_batoms_startup()
print('Successfully setting preference.')
"""

BATOMSPY_TEST = """#!/usr/bin/env batomspy
import bpy
from batoms import Batoms
bpy.ops.batoms.molecule_add()
ch4 = Batoms(label='CH4')
"""

################################################################################
# Section 2.1: shell script (unix only)
################################################################################


# Similar to blender-softwaregl, invoke blender.origin but fix LD
# This script is intended to replace <blender-dir>/blender
BLENDER_REPLACE_SH = """#!/bin/sh
BF_DIST_BIN=$(dirname "$0")
BF_PROGRAM="blender.origin"

LD_LIBRARY_PATH={dyn_lib_path}:${{LD_LIBRARY_PATH}}

if [ -n "$LD_LIBRARYN32_PATH" ]; then
    LD_LIBRARYN32_PATH={dyn_lib_path}:${{LD_LIBRARYN32_PATH}}
fi
if [ -n "$LD_LIBRARYN64_PATH" ]; then
    LD_LIBRARYN64_PATH={dyn_lib_path}:${{LD_LIBRARYN64_PATH}}
fi
if [ -n "$LD_LIBRARY_PATH_64" ]; then
    LD_LIBRARY_PATH_64={dyn_lib_path}:${{LD_LIBRARY_PATH_64}}
fi


export LD_LIBRARY_PATH LD_LIBRARYN32_PATH LD_LIBRARYN64_PATH LD_LIBRARY_PATH_64 LD_PRELOAD

exec "$BF_DIST_BIN/$BF_PROGRAM" ${{1+"$@"}}
"""

# A blender wrapper to be inserted into the conda bin dir
# This is intended to be used for any users, but could be more useful if the user
# want to quickly use `blender` alias inside the conda env
BLENDER_ALIAS_SH = """#!/bin/sh
# This is an alias for blender but with LD_LIBRARY_PATH set
# it is only intended to be used inside the blender conda environment
export LD_LIBRARY_PATH={conda_prefix}/lib:${{LD_LIBRARY_PATH}}
os_name=$(uname -s)

# Add DYLD_LIBRARY_PATH which may be relevant for Dawrin 
# see https://docs.conda.io/projects/conda-build/en/stable/resources/use-shared-libraries.html
if [ "$os_name" = "Darwin" ]
then
    export DYLD_LIBRARY_PATH={conda_prefix}/lib:${{DYLD_LIBRARY_PATH}}
fi

exec "{blender_bin}" ${{1+"$@"}}
"""

# Base for batomspy shabang support
BATOMSPY_SH = """#!/bin/bash
# Usage:
# blenderpy script.py [options]
export LD_LIBRARY_PATH={conda_prefix}/lib:${{LD_LIBRARY_PATH}}
os_name=$(uname -s)

if [ "$os_name" = "Darwin" ]
then
    export DYLD_LIBRARY_PATH={conda_prefix}/lib:${{DYLD_LIBRARY_PATH}}
fi

if [ "$#" -eq  "0" ]
then
    {blender_bin} -b --python-console
else
    {blender_bin} -b --python-exit-code 1 -P $1 ${{@:2}}
fi
"""

################################################################################
# Section 2.2: windows batch files (TBD)
################################################################################


# The directory to move factory python from conda


################################################################################
# Section 3.1: General Utils
################################################################################


def cprint(content, color=None, **kwargs):
    """Color print wrapper"""
    ansi_color = dict(
        HEADER="\033[95m",
        OKBLUE="\033[94m",
        OKGREEN="\033[92m",
        WARNING="\033[93m",
        FAIL="\033[91m",
        ENDC="\033[0m",
        BOLD="\033[1m",
        UNDERLINE="\033[4m",
    )
    if color is None:
        output = content
    elif color in ansi_color.keys() and color != "ENDC":
        output = ansi_color[color] + content + ansi_color["ENDC"]
    else:
        raise ValueError(
            f"Unknown ANSI color name. Allowed values are {list(ansi_color.keys())}"
        )
    print(output, **kwargs)
    return


def _run_process(commands, shell=False, print_cmd=True, cwd=".", capture_output=False):
    """Wrap around subprocess.run
    Returns the process object
    """
    for i in range(len(commands)):
        if isinstance(commands[i], Path):
            commands[i] = commands[i].as_posix()
    full_cmd = " ".join(commands)
    if print_cmd:
        cprint(" ".join(commands))
    if shell is False:
        proc = subprocess.run(
            commands, shell=shell, cwd=cwd, capture_output=capture_output
        )
    else:
        proc = subprocess.run(
            full_cmd, shell=shell, cwd=cwd, capture_output=capture_output
        )
    if proc.returncode == 0:
        return proc
    else:
        raise RuntimeError((f"Running {full_cmd} returned error code {proc.returncode}. \n"
                            "Below are the stdout message: \n"
                            f"{proc.stdout.decode('utf8')}"
                            "\n"
                            "Below are the stderr message:\n"
                            f"{proc.stderr.decode('utf8')}"
                            ))


def _run_blender_multiline_expr(blender_bin, expr, return_process=True, **kwargs):
    """Use blender's interpreter to run multiline
    python expressions. Use `capture_output` to switch whether the output should be printed to stdout
    """
    blender_bin = str(blender_bin)
    tmp_del = False if _get_os_name() in ["windows"] else True
    # The multiline expression must capture
    if "capture_output" not in kwargs:
        kwargs["capture_output"] = True
    with tempfile.NamedTemporaryFile(suffix=".py", delete=tmp_del) as py_file:
        with open(py_file.name, "w") as fd:
            fd.writelines(expr)
        commands = [
            blender_bin,
            "-b",
            "--python-exit-code",
            "1",
            "--python",
            py_file.name,
        ]
        proc = _run_process(commands, print_cmd=False, **kwargs)
    if return_process:
        return proc
    else:
        return None


def _gitclone(workdir=".", version="main", url=REPO_GIT):
    """Make a git clone to the directory
    version can be a branch name or tag name
    return the plugin source path name
    """
    workdir = Path(expanduser(expandvars(workdir)))
    clone_into = workdir / REPO_NAME
    os.makedirs(clone_into, exist_ok=True)
    commands = [
        "git",
        "clone",
        f"{url}",
        clone_into,
    ]
    _run_process(commands)
    cprint(f"Cloned repo into directory {clone_into.as_posix()}", color="OKGREEN")
    _gitcheckout(workdir=clone_into, version=version)
    return clone_into


def _gitcheckout(workdir=".", version="main"):
    """Check a git directory to a specific version"""
    workdir = Path(expanduser(expandvars(workdir)))
    commands = [
        "git",
        "checkout",
        f"{version}",
    ]
    try:
        _run_process(commands, cwd=workdir)
        cprint(f"Checkout version {version}", color="OKGREEN")
    except Exception as e:
        raise RuntimeError(f"Failed to checkout version {version}") from e


################################################################################
# Section 3.2: Tools by parsing blender's commandline output
################################################################################


def _get_default_locations(os_name):
    """Get system specific default install locations of blender.
    Choose multiple possible
    """
    os_name = os_name.lower()
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
                "Automatic discovery of blender on Linux is not supported."
                " Please provide the full path to the blender installation location."
            )
        )
    matches = []
    for l in default_locations[os_name]:
        for version in ALLOWED_BLENDER_VERSIONS:
            true_location = Path(
                expandvars(expanduser(l.format(version=version)))
            ).absolute()
            if true_location.is_dir():
                matches.append(true_location)
    # Multiple blender installation may occur on macos
    if len(matches) > 1:
        cprint(f"Multiple blender installations exists:", color="HEADER")
        for i, m in enumerate(matches):
            cprint(f"{i}: {m.as_posix()}", color="HEADER")
        choice = int(input("Choose one (default 0): ") or "0")
        match = matches[choice]
    elif len(matches) == 1:
        match = matches[0]
    else:
        match = None
    if match is None:
        raise FileNotFoundError(
            (
                f"Cannot find Blender>=3.4 in default installation locations. "
                "Please specify the full path to the blender installation location."
            )
        )
    return match


def _get_blender_bin(os_name, blender_bundle_root):
    """Find the system-dependent blender binary file.
    blender_bundle_root contains the distribution bundle
    and has the pattern <blender_root>/<version>/
    """
    blender_bundle_root = Path(blender_bundle_root).resolve()
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
        help_url = "https://beautiful-atoms.readthedocs.io/en/latest/getting_started/installing/index.html"
        if (os_name == "linux") and (
            "Program Files" in blender_bundle_root.absolute().as_posix()
        ):
            extra_msg = f"""
            It seems you are running the installation script in WSL but points to Blender
            installed on Windows host. Please rerun the script directly on Windows host,
            see instructions at {help_url}
            """
        else:
            extra_msg = f"Please refer to instructions at {help_url}"
        raise FileNotFoundError(
            f"Cannot find blender binary at {blender_bin.resolve().as_posix()}!\n{extra_msg}"
        )
    return blender_bin


def _get_blender_py(blender_bin):
    """Get the blender executable path from current blender binary"""
    blender_bin = str(blender_bin)
    expr = "import sys; print('Python binary: ', sys.executable)"
    proc = _run_blender_multiline_expr(blender_bin, expr)
    output = proc.stdout.decode("utf8")
    cprint(output)
    pat = r"Python\s+binary\:\s+(.*)$"
    py = next(re.finditer(pat, output, re.MULTILINE))[1]
    py_path = Path(py.strip())
    return py_path


def _get_factory_versions(blender_bin):
    """Get the blender version, bundled python and numpy versions
    This is only to be run BEFORE symlinking
    """
    blender_bin = str(blender_bin)
    proc = _run_blender_multiline_expr(blender_bin, BLENDER_CHK_VERSION)
    output = proc.stdout.decode("utf8")
    pat_py_version = r"Python\s+Version\:\s+(\d+\.\d+.\d+)"
    pat_numpy_version = r"Numpy\s+Version\:\s+(\d+\.\d+.\d+)"
    py_version = next(re.finditer(pat_py_version, output))[1]
    numpy_version = next(re.finditer(pat_numpy_version, output))[1]
    return py_version, numpy_version


def _get_blender_version(blender_bin):
    """Parse blender's output to get the X.Y.Z version name"""
    blender_bin = str(blender_bin)
    commands = [blender_bin, "-b"]
    proc = _run_process(commands, shell=False,
                        capture_output=True)
    output = proc.stdout.decode("utf8")
    pat_bl_version = r"Blender\s+(\d+\.\d+\.\d+)"
    bl_version = next(re.finditer(pat_bl_version, output))[1]
    return bl_version


def _blender_enable_plugin(parameters):
    """Use blender's internal libary to enable plugin (and save as user script)"""
    blender_bin = parameters["blender_bin"]
    _run_blender_multiline_expr(blender_bin, BLENDERPY_ENABLE_PLUGIN)
    return


def _blender_disable_plugin(parameters):
    """Use blender's internal libary to disable plugin (and save as user script)"""
    blender_bin = parameters["blender_bin"]
    _run_blender_multiline_expr(blender_bin, BLENDERPY_DISABLE_PLUGIN,)
    print("finish output")
    return


def _blender_test_plugin(parameters):
    if parameters["dependency_only"]:
        cprint("Skip plugin test.", color="WARNING")
        return
    blender_bin = str(parameters["blender_bin"])
    _run_blender_multiline_expr(blender_bin, BLENDERPY_TEST_PLUGIN,)
    return


def _blender_test_uninstall(parameters):
    blender_bin = str(parameters["blender_bin"])
    _run_blender_multiline_expr(blender_bin, BLENDERPY_TEST_UNINSTALL,
                                # capture_output=False
                                )
    return


################################################################################
# Section 3.3: File system tools
################################################################################
def _is_binary_file(filename):
    """Use heurestics to determine if a file is binary
    executable file (containing 0x0) bit
    """
    with open(filename, "rb") as f:
        for i in range(1024):
            try:
                byte = f.read(1)
            except Exception:
                break
            if byte == b"":
                # end of file
                return False
            elif byte == b"\0":
                # null byte found
                return True
    return False


def _is_empty_dir(p):
    """Determin if path has no child dirs"""
    p = Path(p)
    if p.is_dir() is False:
        raise FileNotFoundError(f"{p} is not a valid directory")
    try:
        next(p.rglob("*"))
        stat = False
    except StopIteration:
        stat = True
    return stat


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


def _rename_dir(src, dst):
    """Rename dir from src to dst use os.rename functionality if possible"""
    src = Path(src)
    dst = Path(dst)
    if not src.is_dir():
        raise FileNotFoundError(f"{src} should be an existing directory!")
    if dst.exists():
        if (dst.is_symlink()) or (dst.is_file()):
            os.unlink(dst)
        elif _is_empty_dir(dst):
            os.rmdir(dst)
        else:
            raise OSError(f"Directory {dst.as_posix()} is not empty!")
    try:
        os.rename(src, dst)
    except (OSError, PermissionError) as e:
        raise


def _symlink_dir(src, dst):
    """Make symlink from src to dst.
    If dst is an empty directory or simlink, simply remove and do symlink
    """
    src = Path(src)
    dst = Path(dst)
    if not src.is_dir():
        raise FileNotFoundError(f"{src} should be an existing directory!")
    if dst.exists():
        if (dst.is_symlink()) or (dst.is_file()):
            os.unlink(dst)
        elif _is_empty_dir(dst):
            os.rmdir(dst)
        else:
            raise OSError(f"Directory {dst.as_posix()} is not empty!")
    try:
        os.symlink(src, dst)
    except (OSError, PermissionError) as e:
        raise


################################################################################
# Section 3.4: Conda tools
################################################################################


def _get_conda_variables():
    """Path-objects for conda variables"""
    results = {}
    for key in ("CONDA_PREFIX", "CONDA_PYTHON_EXE", "CONDA_DEFAULT_ENV", "CONDA_EXE"):
        value = os.environ.get(key, "")
        if (value != "") and (key != "CONDA_DEFAULT_ENV"):
            value = Path(value).resolve()
        results[key] = value
    return results


def _is_conda():
    return "CONDA_PREFIX" in os.environ.keys()

def _is_conda_name_abbrev(env_name):
    """Tell if the env_name is a valid conda-env name.
    If not, the env_name may be a full path to the env
    """
    return not any([c in env_name for c in (" ", "/", ":", "#")])


def _find_conda_bin_path(env_name, conda_vars):
    """Return the path of binary search directory ($PATH) of the given conda env"""
    if env_name is None:
        env_name = conda_vars["CONDA_DEFAULT_ENV"]

    # Distinguish if env_name is a prefix or name
    if _is_conda_name_abbrev(env_name):
        commands = [conda_vars["CONDA_EXE"], "run", "-p", env_name, "which", "python"],
    else:
        commands = [conda_vars["CONDA_EXE"], "run", "-n", env_name, "which", "python"],

    
    output = (
        _run_process(
            commands,
            capture_output=True,
        )
        .stdout.decode("utf8")
        .strip()
    )
    bindir = Path(output).parent.resolve()
    return bindir


def _replace_conda_env(python_version=None, numpy_version=None):
    """Replace the env.yml with actual python and numpy versions"""
    blender_py_ver = (
        python_version if python_version is not None else DEFAULT_BLENDER_PY_VER
    )
    blender_numpy_ver = (
        numpy_version if numpy_version is not None else DEFAULT_BLENDER_NUMPY_VER
    )
    env_file_content = ENV_YAML.format(
        blender_py_ver=blender_py_ver, blender_numpy_ver=blender_numpy_ver
    )
    return env_file_content


def _ensure_mamba(conda_vars):
    """Ensure mamba is installed at the base conda environment
    Note: if mamba is not avaivable or install to base is not possible, let user use conda instead
    """
    proc = _run_process(
        [conda_vars["CONDA_EXE"], "list", "-n", "base", "mamba"], capture_output=True
    )
    output = proc.stdout.decode("utf8")
    if "mamba" not in output:
        commands = [
            conda_vars["CONDA_EXE"],
            "install",
            "-y",
            "-n",
            "base",
            "-c",
            "conda-forge",
            "mamba",
        ]
        try:
            _run_process(commands)
        except RuntimeError as e:
            msg = (
                "Failed to install mamba install conda base environment. "
                "You probably don't have write permission. \n"
                "Please consider add --no-mamba to install.py"
            )
            cprint(msg, color="ERROR")
            raise RuntimeError(msg) from e
    # Get the mamba binary in given env
    output = _run_process(
        [conda_vars["CONDA_EXE"].as_posix(), "run", "-n", "base", "which", "mamba"],
        capture_output=True,
    ).stdout.decode("utf8")
    if "ERROR" in output:
        msg = (
            "Cannot find mamba in your conda base environment"
            "Please consider add --no-mamba to install.py"
        )
        raise RuntimeError(output)
    return output.strip()


def _conda_update(
    conda_env_file,
    conda_vars,
    env_name=None,
    python_version=None,
    numpy_version=None,
    backend="mamba",
):
    """Update conda environment using env file.
    If env_name is None, use default env
    If provided python and numpy version, use a temp file to install
    if reinstall_numpy, use pip to reinstall numpy (windows only)
    """
    if env_name is None:
        env_name = conda_vars["CONDA_DEFAULT_ENV"]

    if env_name in ["base"]:
        cprint(
            (
                "Seems you're installing into the base environment. "
                "Installing batoms dependencies may interrupt your base environment. "
            ),
            color="WARNING",
        )
        choice = str(input("Continue? [y/N]") or "N").lower().startswith("y")
        if not choice:
            # TODO: update installation instruction
            cprint(
                "Abort. Please check the installation manual about how to activate an additional conda environment.",
                color="FAIL",
            )
            sys.exit(0)

    # Install from the env.yaml
    cprint("Updating conda environment")
    if backend == "mamba":
        conda_bin = _ensure_mamba(conda_vars)
    else:
        conda_bin = conda_vars["CONDA_EXE"]

    if _is_conda_name_abbrev(env_name):
        commands_prefix = [conda_bin, "env", "update", "-n", env_name]
    else:
        commands_prefix = [conda_bin, "env", "update", "-p", env_name]

    # NamedTemporaryFile can only work on Windows if delete=False
    # see https://stackoverflow.com/questions/55081022/python-tempfile-with-a-context-manager-on-windows-10-leads-to-permissionerror
    tmp_del = False if _get_os_name() in ["windows"] else True
    with tempfile.NamedTemporaryFile(suffix=".yml", delete=tmp_del) as ftemp:
        tmp_yml = ftemp.name
        with open(tmp_yml, "w") as fd:
            env_file_content = _replace_conda_env(python_version, numpy_version)
            fd.writelines(env_file_content.lstrip())

        commands = commands_prefix + [
            "--file",
            tmp_yml,
        ]
        _run_process(commands)
        cprint("Finished install conda packages.", color="OKGREEN")

    return


def _conda_cache_move(condition, conda_vars, blender_python_root):
    """The _conda_cache_move function is DEPRECATED since batoms 2.2.0 as spglib >= 2.0 now has full wheel support on Windows. (Assume no body is using ARM windows 11?)

    Install relevant package (spglib) in conda environment and move to python's site-packages
    This is A DANGEROUS WORKAROUND as it can silently break many things.
    Only to use until the conda environment bug in blender fixed.
    blender_python_root looks like <root>/<3.x>/python
    """
    # Step 1: search latest spglib available for the py version
    conda_bin = conda_vars["CONDA_EXE"]
    cprint(f"Conda bin at: {conda_bin}", color="HEADER")
    commands = [conda_bin, "search", "-c", "conda-forge", str(condition)]
    proc = _run_process(commands, capture_output=True)
    lines = [l for l in proc.stdout.decode("utf-8").split("\n") if len(l) > 1]

    # Choose the newest version
    name, version, build, channel = lines[-1].strip().split()
    conda_url = (
        "https://anaconda.org/conda-forge/{name}/{version}/"
        "download/win-64/{name}-{version}-{build}.tar.bz2"
    ).format(name=name, version=version, build=build)
    # Step 2: do a temp install of spglib into conda environment
    commands = [conda_bin, "install", "--no-deps", conda_url]
    _run_process(commands)
    cprint(f"Installed conda package from {conda_url}", color="OKGREEN")

    # Step 3: copy the site-packages contents
    lib_conda = Path(conda_vars["CONDA_PREFIX"]) / "Lib" / "site-packages"
    lib_pip = Path(blender_python_root) / "lib" / "site-packages"
    match_dirs = lib_conda.glob("spglib*")
    for dir in match_dirs:
        name = dir.name
        shutil.copytree(dir, lib_pip / name, dirs_exist_ok=True)
    return


def _pip_install(blender_py):
    """Temporary workaround for installation on windows and blender>=3.1.0
    Try to install as many components as possible. Need specific version tweaks
    Installation order:
    1. factory numpy -- pinned
    2. ase / scipy / matplotlib / scikit-image / spglib >= 2.0 / pymatgen (>=2022.02) all come with wheel
    3. install openbabel first, if no compiler found, skip
    """
    blender_py = str(blender_py)
    pip_prefix = [blender_py, "-m", "pip"]
    # Step 1: check numpy installation
    commands = pip_prefix + ["show", "numpy"]
    proc = _run_process(commands, shell=False, capture_output=True)
    if "Name: numpy" not in proc.stdout.decode("utf8"):
        raise RuntimeError(
            "Cannot find package numpy. Your bundle python may be corrupt."
        )
    commands = pip_prefix + [
        "install",
        "--no-input",
        "ase>=3.21.0",
        "pymatgen<=2023.3.10",
        "scikit-image",
        "spglib>=2.0.0",
    ]
    proc = _run_process(commands)

    # Step 4: install openbabel (only if compiler exists)
    commands = pip_prefix + ["install", "--no-input", "openbabel"]
    try:
        proc = _run_process(commands)
    except RuntimeError as e:
        cprint(
            (
                "Cannot install openbabel. You need to have a working compiler on windows. "
                "The installation will continue but some functionalities in batoms may not be working."
            ),
            color="WARNING",
        )

    return


def _pip_uninstall(blender_py, conda_vars):
    """uninstall pip components (windows only)"""
    blender_py = str(blender_py)
    pip_prefix = [blender_py, "-m", "pip"]
    # not exhaustive, but should be the most dependent ones
    # TODO: cleanup the dependencies
    commands = pip_prefix + [
        "uninstall",
        "-y",
        "ase",
        "scipy",
        "matplotlib",
        "spglib",
        "scikit-image",
        "plotly",
        "Pillow",
        "openbabel",
        "mpmath",
        "mp_api",
        "monty",
        "msgpack",
        "latexcodec",
        "pybtex",
        "networkx",
        "pandas",
    ]
    _run_process(commands)
    return


def _find_conda_bin_path(env_name, conda_vars):
    """Return the path of binary search directory ($PATH) of the given conda env"""
    if env_name is None:
        env_name = conda_vars["CONDA_DEFAULT_ENV"]

    output = (
        _run_process(
            [
                conda_vars["CONDA_EXE"].as_posix(),
                "run",
                "-n",
                env_name,
                "which",
                "python",
            ],
            capture_output=True,
        )
        .stdout.decode("utf8")
        .strip()
    )
    bindir = Path(output).parent.resolve()
    return bindir


################################################################################
# Section 4.1: operations on blender file structure
################################################################################


def _restore_factory_python(parameters):
    """Restore the factory python directory and remove symlinks
    Logic:
    1)
    """
    # source: <root>/python
    # target: <root>/_old_python
    factory_python_source = parameters["factory_python_source"]
    factory_python_target = parameters["factory_python_target"]

    if factory_python_target.exists():
        # _rename_dir will remove directory if it's empty or a symlink
        try:
            _rename_dir(factory_python_target, factory_python_source)
        except Exception as e:
            raise RuntimeError("Cannot move restore factory python") from e
    # Test if the restored factory python is corrupt?
    if _is_empty_dir(factory_python_source):
        raise RuntimeError("Factory python folder is empty!")
    if factory_python_source.is_symlink():
        raise RuntimeError("Factory python folder is a symlink!")
    try:
        old_py = next(factory_python_source.glob("bin/python*"))
    except StopIteration as e:
        raise RuntimeError("Factory python binary not found!") from e
    # Finally test if the python binary still works
    try:
        _run_process([str(old_py), "-V"])
    except RuntimeError as e:
        raise RuntimeError(
            f"Found factory python at {old_py.as_posix()} "
            "but it's not working. Your installation is corrupted."
        ) from e
    cprint(
        f"Restored {factory_python_target.as_posix()} to {factory_python_source.as_posix()}",
        color="OKGREEN",
    )
    return


def _move_factory_python(parameters):
    """Reverse of _restore_factory_python.
    No sanity check of current python directory
    """
    # blender_root = Path(parameters["blender_root"])
    factory_python_source = parameters["factory_python_source"]
    factory_python_target = parameters["factory_python_target"]

    try:
        _rename_dir(factory_python_source, factory_python_target)
    except Exception as e:
        raise RuntimeError("Cannot move backup factory python") from e
    cprint(
        f"Renamed {factory_python_source.as_posix()} to {factory_python_target.as_posix()}",
        color="OKGREEN",
    )
    return


def _link_conda_env(parameters):
    """Make a symlink from conda environment --> <root>/python"""
    # blender_root = Path(parameters["blender_root"])
    factory_python_source = parameters["factory_python_source"]
    factory_python_target = parameters["factory_python_target"]
    conda_vars = parameters["conda_vars"]
    conda_prefix = conda_vars["CONDA_PREFIX"]
    try:
        _symlink_dir(conda_prefix, factory_python_source)
    except OSError:
        raise
    cprint(
        f"Created symlink {conda_prefix} --> {factory_python_source.as_posix()}",
        color="OKGREEN",
    )
    return




def _install_plugin(parameters):
    """Install plugin into addons_contrib folder
    If in the develop mode, only symlink but not copy
    """
    plugin_path_source = parameters["plugin_source"]
    plugin_path_target = parameters["plugin_target"]
    if parameters["dependency_only"]:
        cprint(
            "Install dependency only. Batoms plugin will not be copied.",
            color="OKGREEN",
        )
        return
    if plugin_path_target.is_dir():
        if not _is_empty_dir(plugin_path_target):
            cprint(
                f"Target plugin installtion directory {plugin_path_target.as_posix()} is not empty.",
                color="WARNING",
            )
            choice = str(input("Overwrite? [y/N]") or "N").lower().startswith("y")
            if not choice:
                cprint("Abort.", color="FAIL")
                sys.exit(0)
        if plugin_path_target.is_symlink():
            os.unlink(plugin_path_target)
        else:
            shutil.rmtree(plugin_path_target)

    # Cases where plugin version should be specified / download from web
    local_plugin_version = parameters["plugin_version"]
    # if the local_plugin_version is not None, download the latest to tmpdir
    if _is_empty_dir(plugin_path_source) or (not plugin_path_source.is_dir()):
        cprint(
            f"Local repo path {plugin_path_source.as_posix()} does not exist. Run git clone.",
        )
        local_plugin_version = "main"
    if local_plugin_version:
        tempdir = tempfile.mkdtemp()
        plugin_path_source = (
            _gitclone(tempdir, local_plugin_version) / DEFAULT_PLUGIN_NAME
        ).resolve()
        # develop option should be disabled
        parameters["develop"] = False

    if parameters["develop"]:
        _symlink_dir(plugin_path_source, plugin_path_target)
        cprint("Installation in development mode!", color="HEADER")
        cprint(
            f"Created symlink {plugin_path_source.as_posix()} --> {plugin_path_target.as_posix()}.",
            color="HEADER",
        )
    else:
        shutil.copytree(plugin_path_source, plugin_path_target)
        cprint(f"Plugin copied to {plugin_path_target.as_posix()}.", color="OKGREEN")
    if "tempdir" in locals():
        shutil.rmtree(tempdir)

    return

def _uninstall_plugin(parameters):
    blender_root = Path(parameters["blender_root"])
    plugin_path_target = blender_root / DEFAULT_PLUGIN_PATH
    if plugin_path_target.is_symlink():
        os.unlink(plugin_path_target)
    elif plugin_path_target.is_dir():
        shutil.rmtree(plugin_path_target)
    else:
        cprint(
            f"Plugin directory {plugin_path_target.as_posix()} does not exist. Ignore.",
            color="WARNING",
        )
    cprint("beautiful-atoms has been removed from blender's addons_contrib folder.", color="OKGREEN")
    return


def _replace_blender_binary(parameters):
    """Replace factory blender binary with a shell script
    that handles LD_LIBRARY_PATH. Only used for linux
    """
    print(parameters)
    blender_bin = parameters["blender_bin"]
    conda_vars = parameters["conda_vars"]
    dyn_lib_path = parameters["blender_root"].resolve() / "python" / "lib"
    if parameters["os_name"] == "linux":
        backup_blender_bin = blender_bin.with_name(BLENDER_BACKUP_PATH)
        if _is_binary_file(blender_bin):
            shutil.move(blender_bin, backup_blender_bin)
            cprint(
                f"Move the blender binary to location {backup_blender_bin}",
                color="HEADER",
            )
        script_path = blender_bin
        with open(script_path, "w") as fd:
            content = BLENDER_REPLACE_SH.format(
                blender_bin=backup_blender_bin.as_posix(),
                conda_prefix=conda_vars["CONDA_PREFIX"].as_posix(),
                dyn_lib_path=dyn_lib_path.as_posix(),
            )
            fd.write(content)
        os.chmod(script_path, 0o755)
        cprint(f"blender script wrapper written to {script_path}", color="OKGREEN")
    else:
        cprint(
            "blender script currently ignored in non-linux platforms", color="HEADER"
        )
    return


def _restore_blender_binary(parameters):
    """Remove the script alias of blender binary"""
    if parameters["os_name"] != "linux":
        return
    blender_bin = parameters["blender_bin"]
    backup_blender_bin = blender_bin.with_name(BLENDER_BACKUP_PATH)
    if backup_blender_bin.exists():
        if not _is_binary_file(backup_blender_bin):
            raise RuntimeError(f"{backup_blender_bin} exists but is not a binary file!")
        if _is_binary_file(blender_bin):
            raise RuntimeError(
                f"Both {blender_bin} and {backup_blender_bin} exist. There is something wrong with your installation!"
            )
        os.remove(blender_bin)
        shutil.move(backup_blender_bin, blender_bin)
        cprint(
            f"Restored the blender binary location {backup_blender_bin}", color="HEADER"
        )
    if not _is_binary_file(blender_bin):
        raise RuntimeError(f"{blender_bin} is corrupted!")
    return


def _create_blender_alias(parameters):
    """Install the alias `blender` script in current conda env."""
    if parameters["os_name"] == "windows":
        cprint("blender alias currently ignored in windows", color="HEADER")
        return
    blender_bin = parameters["blender_bin"]
    conda_vars = parameters["conda_vars"]
    bindir = _find_conda_bin_path(
        env_name=parameters["custom_conda_env"], conda_vars=conda_vars
    )
    script_path = bindir / "blender"
    with open(script_path, "w") as fd:
        content = BLENDER_ALIAS_SH.format(
            blender_bin=blender_bin.as_posix(),
            conda_prefix=conda_vars["CONDA_PREFIX"].as_posix(),
        )
        fd.write(content)
    os.chmod(script_path, 0o755)
    cprint(f"Alias for blender written to {script_path}", color="OKGREEN")

    # Make a test
    
    # with tempfile.NamedTemporaryFile(suffix=".py", delete=False) as tmp_py:
    #     tmp_py = Path(tmp_py.name)
    #     with open(tmp_py, "w") as fd:
    #         fd.write(BATOMSPY_TEST)
    #         os.chmod(tmp_py, 0o755)
    #     _run_process([tmp_py], capture_output=True)
    #     cprint(f"Shebang support for batomspy is now activated",
    #            color="OKGREEN")
    return


def _remove_blender_alias(parameters):
    if parameters["os_name"] == "windows":
        return
    conda_vars = parameters["conda_vars"]
    bindir = _find_conda_bin_path(
        env_name=parameters["custom_conda_env"], conda_vars=conda_vars
    )
    script_path = bindir / "blender"
    if script_path.is_file():
        os.remove(script_path)
        cprint(f"Removed blender alias from {script_path}", color="OKGREEN")
    return


def _create_batomspy(parameters):
    """Install the entrypoint `batomspy` into conda environment's bindir"""
    if parameters["os_name"] == "windows":
        cprint("batoms script currently ignored in windows", color="HEADER")
        return
    blender_bin = parameters["blender_bin"]
    conda_vars = parameters["conda_vars"]
    bindir = _find_conda_bin_path(
        env_name=parameters["custom_conda_env"], conda_vars=conda_vars
    )
    script_path = bindir / "batomspy"
    with open(script_path, "w") as fd:
        content = BATOMSPY_SH.format(
            blender_bin=blender_bin.as_posix(),
            conda_prefix=conda_vars["CONDA_PREFIX"].as_posix(),
        )
        fd.write(content)
    os.chmod(script_path, 0o755)
    cprint(f"batomspy script written to {script_path}", color="OKGREEN")

    # Make a test
    with tempfile.NamedTemporaryFile(suffix=".py", delete=False) as tmp_py:
        tmp_py = Path(tmp_py.name)
        with open(tmp_py, "w") as fd:
            fd.write(BATOMSPY_TEST)
        os.chmod(tmp_py, 0o755)
    _run_process([tmp_py], capture_output=True)
    cprint(f"Shebang support for batomspy is now activated",
           color="OKGREEN")
    return


def _remove_batomspy(parameters):
    if parameters["os_name"] == "windows":
        return
    conda_vars = parameters["conda_vars"]
    bindir = _find_conda_bin_path(
        env_name=parameters["custom_conda_env"], conda_vars=conda_vars
    )
    script_path = bindir / "batomspy"
    if script_path.is_file():
        os.remove(script_path)
        cprint(f"Removed batomspy script from {script_path}", color="OKGREEN")
    return


def _install_dependencies(parameters):
    """Install dependencies from conda or pip"""
    blender_bin = parameters["blender_bin"]
    blender_version = _get_blender_version(blender_bin)
    blender_py = _get_blender_py(blender_bin)
    factory_py_ver, factory_numpy_ver = _get_factory_versions(blender_bin)
    # conda_vars = parameters["conda_vars"]

    if parameters["generate_env_file"] is not None:
        env_file_name = Path(parameters["generate_env_file"])
        with open(env_file_name, "w") as fd:
            fd.writelines(_replace_conda_env(factory_py_ver, factory_numpy_ver))
        cprint(
            f"Conda env exported to {env_file_name.as_posix()}. Exit.", color="OKGREEN"
        )
        sys.exit(0)
    if parameters["use_pip"]:
        _pip_install(blender_py)
        return
    # Rest of the installation will use conda
    if parameters["no_mamba"]:
        backend = "conda"
    else:
        backend = "mamba"
    _conda_update(
        parameters["conda_env_file"],
        parameters["conda_vars"],
        python_version=factory_py_ver,
        numpy_version=factory_numpy_ver,
        backend=backend,
    )
    return

def _uninstall_dependencies(parameters):
       # _old_python not found, ignore
    if not parameters["use_pip"]:
        cprint("Dependencies in conda environment will not be removed.", color="HEADER")
        return

    if factory_python_target.exists():
        raise RuntimeError(
            "Uninstall via pip cannot be performed when bundled python moved to another location. Abort"
        )
        
    # 
    blender_bin = parameters["blender"]
    blender_py = _get_blender_py(blender_bin)
    _pip_uninstall(blender_py, conda_vars)
    cprint("Pip dependencies for beautiful-atoms have been removed", color="OKGREEN")
    return
    
    

def _setup_startup_and_preferences(parameters):
    """Setup blender startup and preferences
    """
    def _blender_set_startup(blender_bin):
        """Use set default startup of batoms"""
        _run_blender_multiline_expr(blender_bin, BLENDERPY_SETTING_STARTUP)
        return


    def _blender_set_preferences(blender_bin):
        """Use set default preferences of batoms"""
        _run_blender_multiline_expr(blender_bin, BLENDERPY_SETTING_PREFERENCES)
        return

    if parameters["use_startup"]:
        _blender_set_startup(blender_bin)
    if parameters["use_preferences"]:
        _blender_set_preferences(blender_bin)
    return


def _print_success_install(parameters):
    """Print message after installation
    """
    print(parameters)
    conda_vars = parameters["conda_vars"]
    env_pref = conda_vars["CONDA_PREFIX"]
    blender_bin = parameters["blender_bin"]
    base_msg = ("beautiful-atoms and its dependencies have been successfully installed in your Blender distribution.\n"
            "If you want to add additional python packages, activate the conda environment first:\n"
            "\n"
            f"$ conda activate --prefix {env_pref.as_posix()}\n"
            "\n"
            "And use conda or pip for the installation.")
    blender_alias_msg = ("\n\n"
                         f"You can now use the alias `blender` for {blender_bin.as_posix()} when the conda environment is activated."
                         )
    batomspy_msg = ("\n\n"
                    f"You can now use `batomspy` in the shebang line of python scripts. For more details please check beautful-atoms' documentation.")
    if parameters["os_name"] != "windows":
        msg = base_msg + blender_alias_msg + batomspy_msg
    else:
        msg = base_msg
    msg += "\n\nHappy coding!"
    cprint(msg, color="OKGREEN")
    return

def _print_success_uninstall(parameters):
    """Print message after installation
    """
    print(parameters)
    conda_vars = parameters["conda_vars"]
    env_pref = conda_vars["CONDA_PREFIX"]
    blender_bin = parameters["blender_bin"]
    base_msg = ("Beautiful-atoms uninstallation finished!")
    pip_msg = ("All dependencies installed by pip have been removed."
               )
    conda_msg = ("Dependencies installed in the conda environment are not removed.\n"
                 "If you don't need to reuse the conda environment, remove it by:\n"
                 "\n"
                 "$ conda deactivate\n"
                 f"$ conda env remove --prefix {env_pref}\n"
                 "\n"
                 )
    if parameters["use_pip"]:
        msg = base_msg + pip_msg
    else:
        msg = base_msg + conda_msg

    cprint(msg, color="OKGREEN")
    return

def _install_version_check(parameters):
    """Pre-installation check of Blender version
    to make sure batoms is compatible
    """
    true_blender_bin = parameters["blender_bin"]
    plugin_version = parameters["plugin_version"]
    # Perform a check on blender_version
    blender_version = _get_blender_version(true_blender_bin)
    if Version(blender_version) >= Version("3.4"):
        cprint("Blender version compatibility OK.", color="OKGREEN")
        return
    
    if plugin_version is None:
        cprint(
            (
                "Warning: support of beautiful-atoms in Blender < 3.4 is deprecated! \n"
                "I will pin the source code of beautiful-atoms to version blender-3.2 branch (793d6f). \n"
                "You can also manually set --plugin-version blender-3.2 to install.py. \n"
                "To use latest features please install Blender >= 3.4. "
            ),
            color="WARNING",
        )
        # This should be the only place parameter changes in-place
        parameters["plugin_version"] = "blender-3.2"
    else:
        cprint(
            "You have specified the plugin_version option, so I suppose you know what you're doing!",
            color="WARNING",
        )
    return
 
def _sanitize_parameters(parameters):
    """Clean up parameters for the installer / uninstaller"""
    parameters["blender_root"] = Path(parameters["blender_root"])
    parameters["blender_bin"] = Path(parameters["blender_bin"])
    parameters["repo_path"] = Path(parameters["repo_path"])
    parameters["conda_env_file"] = parameters["repo_path"] / "env.yml"
    parameters["plugin_source"] = parameters["repo_path"] / DEFAULT_PLUGIN_NAME
    parameters["plugin_target"] = parameters["blender_root"] / DEFAULT_PLUGIN_PATH
    parameters["factory_python_source"] = parameters["blender_root"] / PY_PATH
    parameters["factory_python_target"] = parameters["blender_root"] / PY_BACKUP_PATH
    parameters["conda_vars"] = _get_conda_variables()

    # Only allow --use-pip for now.
    # TODO: check later if conda-env issue is resolved
    if parameters["os_name"] in ["windows"]:
        prev_pip_state = parameters["use_pip"]
        parameters["use_pip"] = True
        if prev_pip_state is not parameters["use_pip"]:
            cprint("For windows installations, I have added --use-pip option for you.", color="HEADER")
        
    # Additional checks
    # 1. are we inside conda?
    if not _is_conda():
        if (parameters["use_pip"] is False):
            cprint(
                "The installation script should be run inside a conda environment or specify the --use-pip option. Abort.",
                color="FAIL",
            )
            sys.exit(0)
            return
        else:
            cprint(("You have specified --use-pip option and conda environment is activated. "
            "Dependencies will be installed directly to blender's python distribution. \n"
            "Some packages like openbabel may require additional setup to install."
                    ), color="WARNING")
    

    # Use `which` to determine if the python executable may
    # conflict with the environment
    if parameters["os_name"] not in ["windows"]:
        check_python_conflict()
    return
    


def install(parameters):
    """Link current conda environment to blender's python root
    Copy batoms plugin under repo_path to plugin directory
    """
    parameters = parameters.copy()
    _sanitize_parameters(parameters)
    _install_version_check(parameters)
    # Linux-only. Replace blender binary with  LD_LIBRARY_PATH awareness
    _restore_blender_binary(parameters)
    # TODO: make paths in binary script relative
    _replace_blender_binary(parameters)
    _restore_factory_python(parameters)
    _move_factory_python(parameters)
    _link_conda_env(parameters)
    # TODO: condition about generate env file?
    # TODO: finish the conda part in between
    _install_dependencies(parameters)
    _install_plugin(parameters)
    _setup_startup_and_preferences(parameters)
    _blender_enable_plugin(parameters)
    _blender_test_plugin(parameters)
    _create_blender_alias(parameters)
    _create_batomspy(parameters)
    _print_success_install(parameters)
    return


def uninstall(parameters):
    """Remove the plugin from target_directory and restore"""
    parameters = parameters.copy()
    _sanitize_parameters(parameters)
    _blender_disable_plugin(parameters)
    _uninstall_plugin(parameters)
    _uninstall_dependencies(parameters)
    _remove_blender_alias(parameters)
    _remove_batomspy(parameters)
    _restore_factory_python(parameters)
    _restore_blender_binary(parameters)
    _blender_test_uninstall(parameters)
    _print_success_uninstall(parameters)
    return



def check_python_conflict():
    """Detect if the python interpreter used to run the installation script
    is the same as the python executable in current environment
    """
    current_py = sys.executable
    env_py = shutil.which("python")
    if env_py is None:
        # Does not have a python interpreter in current environment, safe.
        return
    if Path(current_py).resolve().samefile(Path(env_py)):
        cprint(
            (
                "You're running install.py script using python interpreter:\n"
                f"{current_py}\n"
                "It may be updated during the install process and causing issues. "
                "We recommend using another python interpreter for install.py, such as: \n"
                "$CONDA_PYTHON_EXE install.py [options]"
            ),
            color="WARNING",
        )
        choice = str(input("Continue? [y/N]") or "N").lower().startswith("y")
        if not choice:
            # TODO: update installation instruction
            cprint("Abort.", color="FAIL")
            sys.exit(0)
        else:
            return
    else:
        return





def main():
    import argparse

    curdir = Path(__file__).parent.resolve()

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
        "-t",
        "--plugin-version",
        default=None,
        help=(
            "Plugin version (in the form 'vX.Y.Z') or git hash tag. "
            "Note this option, if not None, discards the values from "
            "--local-repo-path and --develop."
        ),
    )
    parser.add_argument(
        "-p",
        "--local-repo-path",
        default=curdir,
        help=(
            "Path to root of the beautiful-atoms repo. "
            "Default is the same level as install.py. "
            "If using install.py alone, download the latest version from github and install."
        ),
    )
    parser.add_argument(
        "--uninstall", action="store_true", help="Uninstall plugin in blender_root"
    )
    parser.add_argument(
        "--use-startup",
        action="store_true",
        default=False,
        help="Use the default startup.",
    )
    parser.add_argument(
        "--use-preferences",
        action="store_true",
        default=False,
        help="Use the default preferences.",
    )
    parser.add_argument(
        "-n",
        "--conda-env-name",
        default=None,
        help="Conda environment to install dependencies other than current environment.",
    )
    parser.add_argument(
        "--use-pip",
        action="store_true",
        help="Use pip install instead of conda environment. Only recommended on windows.",
    )
    parser.add_argument(
        "--dependency-only",
        action="store_true",
        help=(
            "Install dependencies but not copying batoms plugin. "
            "Should only be used for building docker images."
        ),
    )
    parser.add_argument(
        "--generate-env-file",
        nargs="?",
        const="env.yml",
        help=(
            "Only print the dependency as a env.yml to local dir without any installation."
        ),
    )
    parser.add_argument(
        "--develop",
        action="store_true",
        help=(
            "Development mode. Symlink the local batoms directory to plugin folder in blender."
            "After such installation, local change of batoms source code will have immediate effect "
            "without running install.py again."
        ),
    )
    parser.add_argument(
        "--no-mamba",
        action="store_true",
        help=("Use the default conda backend instead of mamba for version resolver"),
    )
    parser.add_argument(
        "--no-blender-alias",
        action="store_true",
        help=("Do not install the blender alias in conda env's root. Only relevant for unix systems."),
    )
    args = parser.parse_args()
    os_name = _get_os_name()

    
    if (os_name not in ["windows"]) and (args.uninstall is False):
        check_python_conflict()

    # Use blender_root supplied to commandline if exists, otherwise try to guess
    if args.blender_root:
        true_blender_root = Path(expanduser(expandvars(args.blender_root)))
    else:
        true_blender_root = _get_default_locations(os_name)
    true_blender_bin = _get_blender_bin(os_name, true_blender_root)
    cprint(f"Found blender binary at {true_blender_bin.as_posix()}", color="OKGREEN")
    cprint(
        f"      blender bundle root at {true_blender_root.as_posix()}", color="OKGREEN"
    )

    

    # Parameters can be provided to install / uninstall methods at this time
    parameters = dict(
        blender_root=true_blender_root,
        blender_bin=true_blender_bin,
        os_name=os_name,
        use_pip=args.use_pip,
        repo_path=Path(expanduser(expandvars(args.local_repo_path))),
        custom_conda_env=args.conda_env_name,
        dependency_only=args.dependency_only,
        generate_env_file=args.generate_env_file,
        use_startup=args.use_startup,
        use_preferences=args.use_preferences,
        develop=args.develop,
        no_mamba=args.no_mamba,
        no_blender_alias=args.no_blender_alias,
        plugin_version=args.plugin_version,
    )
    if not args.uninstall:
        install(parameters)
    else:
        uninstall(parameters)


if __name__ == "__main__":
    main()
