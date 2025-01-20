### Beautiful Atoms
[![Build beautiful_atoms main image](https://github.com/beautiful-atoms/beautiful-atoms/actions/workflows/build_main_image.yml/badge.svg)](https://github.com/beautiful-atoms/beautiful-atoms/actions/workflows/build_main_image.yml)
[![Test batoms blender plugin](https://github.com/beautiful-atoms/beautiful-atoms/actions/workflows/batoms_plugin_test.yaml/badge.svg)](https://github.com/beautiful-atoms/beautiful-atoms/actions/workflows/batoms_plugin_test.yaml)

Batoms is a Python package for editing and rendering atoms and molecule objects using Blender. A Python interface that allows for automating workflows.

Features:

* Model: space-filling, ball-stick, polyhedral, cavity and so on.
* Supported File type: cif, xyz, cube, pdb, json, VASP-out and so on.
* Supported structure: ASE and Pymatgen
* Volumetric data (Isosurface)
* Ribbon diagram for protein
* Site occupancy
* Animation
* GUI
* Support periodic boundary conditions
* Support fetch structures from MaterialProject, Pubchem and RSCB
* ``Flexible``: Python script, run interactively or in the background.
* ``High quality rendering``:  3D models
* ``Free, Open Source``: Easy to download and install.
* ``Cross-platform``: (Linux, Windows, macOS)


### How to use

Please visit: https://beautiful-atoms.readthedocs.io/en/latest/

### View and edit structure in Jupyter Notebook
Another package, [weas-widget](https://github.com/superstar54/weas-widget), allows users to visualize and edit with atomistic structures in Jupyter Notebook.

### How to contribute

#### Editor
We recommend using [Visual Studio Code](https://code.visualstudio.com/) with the [Blender extension](https://github.com/JacquesLucke/blender_vscode).

#### Test
We recommend using [pytest-blender](https://pypi.org/project/pytest-blender/). To run the tests, run:

```console
pip install pytest-blender
blender_python="$(pytest-blender)"
$blender_python -m ensurepip
# Install the development dependencies:
$blender_python -m pip install -r test-requirements.txt
cd tests
pytest
```

#### Pre-commit
To install the pre-commit hooks, run:

```console
$ pre-commit install
```


### Extension

```
pip download ase --dest ./batoms/wheels --only-binary=:all: --python-version=3.11 --platform=manylinux_2_17_x86_64
pip download ase --dest ./batoms/wheels --only-binary=:all: --python-version=3.11 --platform=win_amd64
pip download ase --dest ./batoms/wheels --only-binary=:all: --python-version=3.11 --platform=macosx_12_0_arm64
```

```
blender --command extension build --source-dir batoms --output-dir .
```

### Contact
* Xing Wang  <xingwang1991@gmail.com>
