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
from distutils import command
import os
import re
from os.path import expanduser, expandvars
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from distutils.version import LooseVersion
from threading import local
import stat

# TODO: allow version control
# TODO: windows privilege issue
# TODO: complete install tutorial about the env variables
# TODO: unit test workflow
# TODO: python && numpy version match
# sys.version --> python version
# numpy.__version__ --> numpy


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
  - spglib
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

BLENDERPY_SETTING_PREFERENCES = f"""
print('start setting preferences')
import bpy
bpy.ops.batoms.use_batoms_preference()
print('Successfully setting and preference.')
"""

BLENDERPY_SETTING_STARTUP = f"""
print('start setting startup')
import bpy
bpy.ops.batoms.use_batoms_startup()
print('Successfully setting preference.')
"""

# wrapper for blender script.
# Replace the value {blender_bin} at installation runtime
BATOMSPY_SH = """#!/bin/bash
# Usage:
# blenderpy script.py [options]
if [ "$#" -eq  "0" ]
then
    {blender_bin} -b --python-console
else
    {blender_bin} -b --python-exit-code 1 -P $1 ${{@:2}}
fi
"""

BATOMSPY_TEST = """#!/usr/bin/env batomspy
import bpy
from batoms import Batoms
bpy.ops.batoms.molecule_add()
ch4 = Batoms(label='CH4')
"""

# The directory to move factory python from conda

repo_name = os.environ.get("GITHUB_REPO", DEFAULT_REPO_NAME)
account_name = os.environ.get("GITHUB_ACCOUNT", DEFAULT_GITHUB_ACCOUNT)

repo_git = f"https://github.com/{account_name}/{repo_name}.git"
install_script_path = (
    f"https://raw.githubusercontent.com/{account_name}/{repo_name}/main/install.py"
)

# Minimal class for terminal ANSI color output


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
            f"Cannot find blender binary at {blender_bin.as_posix()}!\n{extra_msg}"
        )
    return blender_bin


def _get_blender_py(blender_bin):
    """Get the blender executable path from current blender binary"""
    blender_bin = str(blender_bin)
    commands = [
        blender_bin,
        "-b",
        "--python-expr",
        "import sys; print('Python binary: ', sys.executable)",
    ]
    proc = _run_process(commands, shell=False, capture_output=True)
    output = proc.stdout.decode("utf8")
    cprint(output)
    pat = r"Python\s+binary\:\s+(.*)$"
    py = next(re.finditer(pat, output, re.MULTILINE))[1]
    py_path = Path(py.strip())
    return py_path


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


def _get_factory_versions(blender_bin):
    """Get the blender version, bundled python and numpy versions
    This is only to be run BEFORE symlinking
    """
    blender_bin = str(blender_bin)
    # First get blender python version
    commands = [blender_bin, "-b", "--python-expr", BLENDER_CHK_VERSION]
    proc = _run_process(commands, shell=False, capture_output=True)
    output = proc.stdout.decode("utf8")
    pat_py_version = r"Python\s+Version\:\s+(\d+\.\d+.\d+)"
    pat_numpy_version = r"Numpy\s+Version\:\s+(\d+\.\d+.\d+)"
    py_version = next(re.finditer(pat_py_version, output))[1]
    numpy_version = next(re.finditer(pat_numpy_version, output))[1]
    return py_version, numpy_version


def _get_blender_version(blender_bin):
    """Parse blender's output to get the X.Y.Z version name"""
    blender_bin = str(blender_bin)
    # First get blender python version
    commands = [blender_bin, "-b"]
    proc = _run_process(commands, shell=False, capture_output=True)
    output = proc.stdout.decode("utf8")
    pat_bl_version = r"Blender\s+(\d+\.\d+\.\d+)"
    bl_version = next(re.finditer(pat_bl_version, output))[1]
    return bl_version


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


def _run_process(commands, shell=False, print_cmd=True, cwd=".", capture_output=False):
    """Wrap around subprocess.run
    Returns the process object
    """
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
        raise RuntimeError(f"Running {full_cmd} returned error code {proc.returncode}")


def _run_blender_multiline_expr(blender_bin, expr, **kwargs):
    blender_bin = str(blender_bin)
    tmp_del = False if _get_os_name() in ["windows"] else True
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
        _run_process(commands, print_cmd=False, **kwargs)
    return


def _blender_enable_plugin(blender_bin):
    """Use blender's internal libary to enable plugin (and save as user script)"""
    _run_blender_multiline_expr(blender_bin, BLENDERPY_ENABLE_PLUGIN)
    return


def _blender_disable_plugin(blender_bin):
    """Use blender's internal libary to disable plugin (and save as user script)"""
    _run_blender_multiline_expr(blender_bin, BLENDERPY_DISABLE_PLUGIN)
    return


def _blender_test_plugin(parameters):
    if parameters["dependency_only"] is False:
        blender_bin = str(parameters["blender_bin"])
        _run_blender_multiline_expr(blender_bin, BLENDERPY_TEST_PLUGIN)
    else:
        cprint("Skip plugin test.", color="WARNING")
    return


def _blender_test_uninstall(parameters):
    blender_bin = str(parameters["blender_bin"])
    _run_blender_multiline_expr(blender_bin, BLENDERPY_TEST_UNINSTALL)
    return


def _gitclone(workdir=".", version="main", url=repo_git):
    """Make a git clone to the directory
    version can be a branch name or tag name
    return the plugin source path name
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
    cprint(f"Cloned repo into directory {clone_into.as_posix()}", color="OKGREEN")
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


def _rename_dir(src, dst):
    """Rename dir from scr to dst use os.rename functionality if possible"""
    src = Path(src)
    dst = Path(dst)
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
        # TODO: do something here
        raise


def _symlink_dir(src, dst):
    """Make symlink from src to dst.
    If dst is an empty directory or simlink, simply remove and do symlink
    """
    src = Path(src)
    dst = Path(dst)
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
        # TODO: do something here
        raise


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
            msg = ("Failed to install mamba install conda base environment. "
            "You probably don't have write permission. \n"
            "Please consider add --no-mamba to install.py"
            )
            cprint(msg,
            color="ERROR")
            raise RuntimeError(msg) from e
    # Get the mamba binary in given env
    output = _run_process(
        [conda_vars["CONDA_EXE"], "run", "-n", "base", "which", "mamba"],
        capture_output=True,
    ).stdout.decode("utf8")
    if "ERROR" in output:
        msg = ("Cannot find mamba in your conda base environment"
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
        mamba_bin = _ensure_mamba(conda_vars)
        # This is dangerous running outside conda env. Consider check the base first
        commands_prefix = [
            mamba_bin,
            "env",
            "update",
            "-n",
            env_name,
        ]
    else:
        commands_prefix = [conda_vars["CONDA_EXE"], "env", "update", "-n", env_name]

    # NamedTemporaryFile can only work on Windows if delete=False
    # see https://stackoverflow.com/questions/55081022/python-tempfile-with-a-context-manager-on-windows-10-leads-to-permissionerror
    tmp_del = False if _get_os_name() in ["windows"] else True
    with tempfile.NamedTemporaryFile(suffix=".yml", delete=tmp_del) as ftemp:
        tmp_yml = ftemp.name
        with open(tmp_yml, "w") as fd:
            env_file_content = _replace_conda_env(python_version, numpy_version)
            fd.writelines(env_file_content.lstrip())

        commands = commands_prefix + [
            # conda_vars["CONDA_EXE"],
            # "env",
            # "update",
            # "-n",
            # env_name,
            "--file",
            tmp_yml,
            # conda_env_file.as_posix(),
        ]
        _run_process(commands)
        cprint("Finished install conda packages.", color="OKGREEN")

    return


def _conda_cache_move(condition, conda_vars, blender_python_root):
    """Install relevant package (spglib) in conda environment and move to python's site-packages
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


def _pip_install(blender_py, blender_python_root, factory_py_ver, conda_vars):
    """Temporary workaround for installation on windows and blender>=3.1.0
    Try to install as many components as possible. Need specific version tweaks
    Installation order:
    1. factory numpy -- pinned
    2. install spglib with pip first, if no compiler found, try building from conda-forge's distribution
    3. ase / scipy / matplotlib / scikit-image all come with wheel
    4. install pymatgen (>=2022.02)
    5. install openbabel first, if no compiler found, skip
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

    # Step 2: install spglib
    commands = pip_prefix + ["install", "--no-input", "spglib"]
    try:
        proc = _run_process(commands, shell=False, capture_output=True)
    except RuntimeError:
        cprint(
            "Building spglib from source failed. We'll try to install from conda-distruted lib.",
            color="WARNING",
        )
        # abbrevate version, i.e. 3.10 --> py310
        abbrev_py_ver = "py" + "".join(factory_py_ver.split(".")[:2])
        condition = f"spglib=*={abbrev_py_ver}*"
        # cprint(condition)
        _conda_cache_move(
            condition=condition,
            conda_vars=conda_vars,
            blender_python_root=blender_python_root,
        )
    # check spglib installation
    commands = pip_prefix + ["show", "spglib"]
    proc = _run_process(commands, shell=False, capture_output=True)
    if "Name: spglib" not in proc.stdout.decode("utf8"):
        # TODO: improve error msg
        raise RuntimeError("Spglib installation failed.")

    # Step 3: install ase pymatgen etc.
    commands = pip_prefix + [
        "install",
        "--no-input",
        "ase>=3.21.0",
        "pymatgen<=2022.03",
        "scikit-image",
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
        "monty",
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
            [conda_vars["CONDA_EXE"], "run", "-n", env_name, "which", "python"],
            capture_output=True,
        )
        .stdout.decode("utf8")
        .strip()
    )
    bindir = Path(output).parent.resolve()
    return bindir


def install(parameters):
    """Link current conda environment to blender's python root
    Copy batoms plugin under repo_path to plugin directory
    """
    local_parameters = parameters.copy()
    blender_root = Path(local_parameters["blender_root"])
    # cprint("blender_root", blender_root)
    blender_bin = Path(local_parameters["blender_bin"])
    repo_path = Path(local_parameters["repo_path"])
    conda_env_file = repo_path / "env.yml"
    # plugin_path_source: dir of plugin in github repo
    # plugin_path_target: dir of plugin under blender root
    plugin_path_source = (repo_path / DEFAULT_PLUGIN_NAME).resolve()
    plugin_path_target = (blender_root / DEFAULT_PLUGIN_PATH).resolve()

    # factory_python_source: (presumably) python shipped with blender
    # factory_python_target: dir to move factory python
    factory_python_source = blender_root / PY_PATH
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
        # Empty symlink factory python might be mistakenly created
        # just delete them
        exception_corrupt = False
        if _is_empty_dir(factory_python_target) or factory_python_target.is_symlink():
            exception_corrupt = True
        else:
            cprint(f"Target path {factory_python_target.as_posix()} exists.")
            try:
                old_py = next(factory_python_target.glob("bin/python*"))
            except StopIteration:
                cprint("Factory python binary not found.", color="WARNING")
                old_py = None
                exception_corrupt = True
            if old_py:
                try:
                    _run_process([str(old_py), "-V"])
                except RuntimeError:
                    cprint(
                        (
                            f"Found factory python at {old_py.as_posix()} "
                            "but it's not working"
                        ),
                        color="WARNING",
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
            # os.rename(factory_python_target, factory_python_source)
            shutil.move(factory_python_target, factory_python_source)
            cprint(
                f"Renamed {factory_python_target.as_posix()} to {factory_python_source.as_posix()}",
                color="OKGREEN",
            )

    # Step 1.1 check python and numpy version
    (
        factory_py_ver,
        factory_numpy_ver,
    ) = _get_factory_versions(blender_bin)
    blender_version = _get_blender_version(blender_bin)
    blender_py = _get_blender_py(blender_bin)

    # User warning for blender_version < 3.4 should be already resolved outside install(parameter)
    # only print the warning again
    # Compare version
    # if LooseVersion(str(version)) < LooseVersion(str(MIN_BLENDER_VER)):
    # raise ValueError(
    # f"Blender version {version} is not supported. Minimal requirement is {MIN_BLENDER_VER}"
    # )
    if LooseVersion(blender_version) in LooseVersion("3.4"):
        cprint("Warning: support for beautiful-atoms in Blender <3.4 is deprecated!",
        color="WARNING")

    print(blender_version, factory_py_ver, factory_numpy_ver)

    # If need to output the env yaml only, do it now
    # cprint(parameters["generate_env_file"])
    if local_parameters["generate_env_file"] is not None:
        env_file_name = Path(local_parameters["generate_env_file"])
        with open(env_file_name, "w") as fd:
            fd.writelines(_replace_conda_env(factory_py_ver, factory_numpy_ver))
        cprint(
            f"Conda env exported to {env_file_name.as_posix()}. Exit.", color="OKGREEN"
        )
        sys.exit(0)

    # Step 2. cond 1: use pip
    if local_parameters["use_pip"]:
        _pip_install(blender_py, factory_python_source, factory_py_ver, conda_vars)

    else:
        # Step 2: cond2: rename soruce to target
        if factory_python_target.exists():
            raise RuntimeError(
                f"Something wrong. {factory_python_target.as_posix()} still exists."
            )
        if factory_python_source.is_symlink():
            os.unlink(factory_python_source)
            cprint(f"{factory_python_source.as_posix()} is already a symlink and will be unlinked. No backup will be made.", color="WARNING")
        elif not factory_python_source.is_dir():
            cprint(f"{factory_python_source.as_posix()} does not exist. Will not move.", color="WARNING")
        else:
            shutil.move(factory_python_source, factory_python_target)
            cprint(
                f"Renamed {factory_python_source.as_posix()} to {factory_python_target.as_posix()}",
                color="OKGREEN",
            )
        # Step 2-1: link the conda prefix of current environment
        conda_prefix = Path(conda_vars["CONDA_PREFIX"]).resolve()
        # Should not happen but just in case
        if factory_python_source.is_symlink():
            os.unlink(factory_python_source)
        os.symlink(conda_prefix, factory_python_source)
        cprint(
            f"Created symlink {conda_prefix} --> {factory_python_source.as_posix()}",
            color="OKGREEN",
        )

        # Give a warning about conda env
        # TODO: allow direct install into another environment
        if local_parameters["no_mamba"]:
            backend = "conda"
        else:
            backend = "mamba"
        _conda_update(
            conda_env_file,
            conda_vars,
            python_version=factory_py_ver,
            numpy_version=factory_numpy_ver,
            backend=backend,
        )

    # Step 4: install plugin
    if local_parameters["dependency_only"]:
        cprint(
            "Install dependency only. Batoms plugin will not be copied.",
            color="OKGREEN",
        )
        return

    if not _is_empty_dir(plugin_path_target):
        cprint(
            f"Target plugin installtion directory {plugin_path_target.as_posix()} is not empty.",
            color="WARNING",
        )
        choice = str(input("Overwrite? [y/N]") or "N").lower().startswith("y")
        if not choice:
            cprint("Abort.", color="FAIL")
            sys.exit(0)
        else:
            if plugin_path_target.is_symlink():
                os.unlink(plugin_path_target)
            else:
                shutil.rmtree(plugin_path_target)
    # Step 5. Checkout the plugin version
    local_plugin_version = local_parameters["plugin_version"]
    # if the local_plugin_version is not None, download the latest to tmpdir
    if local_plugin_version:
        tempdir = tempfile.mkdtemp()
        plugin_path_source = (
            _gitclone(tempdir, local_plugin_version) / DEFAULT_PLUGIN_NAME
        ).resolve()
        # cprint(plugin_path_source)
        local_parameters["develop"] = False
    # if version is None and plugin_path_source does not exist, download the latest to tmpdir
    elif _is_empty_dir(plugin_path_source) or (not plugin_path_source.is_dir()):
        cprint(
            f"Local repo path {plugin_path_source.as_posix()} does not exist. Run git clone.",
        )
        tempdir = tempfile.mkdtemp()
        local_plugin_version = "main"
        plugin_path_source = (
            _gitclone(tempdir, local_plugin_version) / DEFAULT_PLUGIN_NAME
        ).resolve()
        # cprint(plugin_path_source)
        local_parameters["develop"] = False
    else:
        pass

    # Choose whether to copy the plugin path or symlink it

    if local_parameters["develop"]:
        os.symlink(plugin_path_source, plugin_path_target)
        cprint("Installation in development mode!", color="HEADER")
        cprint(
            f"Created symlink {plugin_path_source.as_posix()} --> {plugin_path_target.as_posix()}.",
            color="HEADER",
        )
    else:
        shutil.copytree(plugin_path_source, plugin_path_target)
        cprint(f"Plugin copied to {plugin_path_target.as_posix()}.", color="OKGREEN")
    _blender_enable_plugin(blender_bin)
    _blender_test_plugin(local_parameters)

    if local_parameters["use_startup"]:
        _blender_set_startup(blender_bin)
    if local_parameters["use_preferences"]:
        _blender_set_preferences(blender_bin)

    # Add the batomspy script (for unix systems only)
    if local_parameters["os_name"] != "windows":
        bindir = _find_conda_bin_path(
            env_name=local_parameters["custom_conda_env"], conda_vars=conda_vars
        )
        script_path = bindir / "batomspy"
        with open(script_path, "w") as fd:
            content = BATOMSPY_SH.format(blender_bin=blender_bin.as_posix())
            fd.write(content)
        os.chmod(script_path, 0o755)
        cprint(f"batomspy script written to {script_path}", color="OKGREEN")
    else:
        cprint("batoms script currently ignored in windows", color="WARNING")

    if local_parameters["os_name"] != "windows":
        with tempfile.NamedTemporaryFile(suffix=".py", delete=False) as tmp_py:
            tmp_py = Path(tmp_py.name)
            with open(tmp_py, "w") as fd:
                fd.write(BATOMSPY_TEST)
            os.chmod(tmp_py, 0o755)
        _run_process([tmp_py.as_posix()], capture_output=True)
        cprint(f"Shebang support for batomspy is now activated", color="OKGREEN")

    cprint(
        (
            "Beautiful-atoms and its dependencies have been successfully installed in your Blender distribution.\n"
            "If you want to add additional python packages, activate the conda environment first:\n"
            "\n"
            "$ conda activate --prefix {env_pref}\n"
            "\n"
            "And use conda or pip for the installation. Happy coding!"
        ).format(env_pref=conda_vars["CONDA_PREFIX"]),
        color="OKGREEN",
    )
    return


def uninstall(parameters):
    """Remove the plugin from target_directory and restore"""
    local_parameters = parameters.copy()
    blender_root = Path(local_parameters["blender_root"])
    blender_bin = Path(local_parameters["blender_bin"])
    plugin_path_target = blender_root / DEFAULT_PLUGIN_PATH
    factory_python_source = blender_root / PY_PATH
    factory_python_target = blender_root / PY_BACKUP_PATH

    conda_vars = _get_conda_variables()

    # Remove the batomspy script
    bindir = _find_conda_bin_path(
        env_name=local_parameters["custom_conda_env"], conda_vars=conda_vars
    )
    script_path = bindir / "batomspy"
    if script_path.is_file():
        os.remove(script_path)
        cprint(f"Removed batomspy script from {script_path}", color="OKGREEN")

    _blender_disable_plugin(blender_bin)
    # if parameters["use_pip"]:
    #     _pip_uninstall()
    if plugin_path_target.is_symlink():
        os.unlink(plugin_path_target)
    elif plugin_path_target.is_dir():
        shutil.rmtree(plugin_path_target)
    else:
        cprint(
            f"Plugin directory {plugin_path_target.as_posix()} does not exist. Ignore.",
            color="WARNING",
        )

    # _old_python not found, ignore
    if not factory_python_target.exists():
        if local_parameters["use_pip"]:
            blender_py = _get_blender_py(blender_bin)
            _pip_uninstall(blender_py, conda_vars)
            cprint(
                "Beautiful-atoms and its dependencies have been successfully uninstalled!",
                color="OKGREEN",
            )
        else:
            cprint(
                f"Backup of factory blender python path {factory_python_target.as_posix()} does not exist. Ignore.",
                color="WARNING",
            )
            cprint(
                "It seems beautiful-atoms has already been uninstalled from your Blender distribution.",
                color="OKGREEN",
            )
        return
    else:
        if local_parameters["use_pip"]:
            raise RuntimeError(
                "Uninstall via pip cannot be performed when bundled python moved to another location. Abort"
            )
        if factory_python_source.is_dir():
            if not _is_empty_dir(factory_python_source):
                if factory_python_source.is_symlink():
                    os.unlink(factory_python_source)
                elif factory_python_source.is_dir():
                    cprint(
                        f"Current blender python path is not a symlink.", color="HEADER"
                    )
                    overwrite = (
                        str(input("Overwrite? [y/N]") or "N").lower().startswith("y")
                    )
                    if overwrite:
                        shutil.rmtree(factory_python_source)
                    else:
                        cprint("Ignored.", color="FAIL")
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
    _blender_test_uninstall(local_parameters)
    cprint(
        (
            "Beautiful-atoms uninstallation finished.\n"
            "If you don't need to reuse the conda environment, remove it by:\n"
            "\n"
            "$ conda deactivate\n"
            "$ conda env remove --prefix {env_pref}\n"
            "\n"
        ).format(env_pref=conda_vars["CONDA_PREFIX"]),
        color="OKGREEN",
    )
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


def _blender_set_startup(blender_bin):
    """Use set default startup of batoms"""
    _run_blender_multiline_expr(blender_bin, BLENDERPY_SETTING_STARTUP)
    return


def _blender_set_preferences(blender_bin):
    """Use set default preferences of batoms"""
    _run_blender_multiline_expr(blender_bin, BLENDERPY_SETTING_PREFERENCES)
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
    # TODO: to be implemented later
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
    args = parser.parse_args()
    cprint(args)
    os_name = _get_os_name()

    # When installing on linux / macos using conda method, makesure python interpreters do not conflict
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

    # Perform a check on blender_version
    blender_version = _get_blender_version(true_blender_bin)
    if LooseVersion(blender_version) < LooseVersion("3.4"):
        # The warning 
        if args.plugin_version is None:
            cprint(("Warning: support of beautiful-atoms in Blender < 3.4 is deprecated! "
            "I will pin the source code of beautiful-atoms to version 793d6f (blender-3.2 branch)."
            "To use latest features please install Blender >= 3.4."),
            color="WARNING"
            )
            plugin_version = "793d6f"
        else:
            cprint("You have specified the plugin_version option, so I suppose you know what you're doing!")
            plugin_version = args.plugin_version
    else:
        plugin_version = args.plugin_version

    # Parameters can be provided to install / uninstall methods at this time
    # Do not process any information regarding blender version / python version
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
        plugin_version=plugin_version,
    )

    # Uninstallation does not need information about current environment
    if args.uninstall:
        if (parameters["os_name"] in ["windows"]) and (parameters["use_pip"] is False):
            cprint(
                (
                    "To uninstall batoms depedencies on windows, please add the --use-pip tag to the script."
                ),
                color="FAIL",
            )
            sys.exit(0)
        uninstall(parameters)
        return

    # Cannot install without conda if pip install is enabled
    if (not is_conda()) and (args.use_pip is False):
        cprint(
            "The installation script should be run inside a conda environment. Abort.",
            color="FAIL",
        )
        sys.exit(0)

    # Sanity check. pip install only recommended for windows
    if args.use_pip:
        if os_name not in ["windows"]:
            cprint(
                (
                    "Install dependencies via pip is only recommended for windows."
                    " Please remove the --use-pip flag and try again"
                ),
                color="FAIL",
            )
            sys.exit(0)
    elif os_name in ["windows"]:
        cprint(
            (
                "Conda install currently not working for Blender>=3.0 "
                "due to a bug in anaconda https://github.com/ContinuumIO/anaconda-issues/issues/11994.\n"
                "Please add the --use-pip flag to installation script. "
            ),
            color="FAIL",
        )
        sys.exit(0)

    install(parameters)


if __name__ == "__main__":
    main()
