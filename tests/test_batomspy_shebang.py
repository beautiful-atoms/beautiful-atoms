import subprocess
import pytest
import os
import shutil
from pathlib import Path
import tempfile

BATOMSPY_TEST = """#!/usr/bin/env batomspy
import bpy
from batoms import Batoms
from ase.build import molecule
ch4 = Batoms(from_ase=molecule('CH4'))
ch4.render.samples = 1
ch4.render.resolution = [10, 10]
ch4.get_image(engine="cycles")
"""


def _run_process(commands, shell=False, print_cmd=True, cwd=".", capture_output=False):
    """Wrap around subprocess.run
    Returns the process object. Copied from install.py
    """
    full_cmd = " ".join(commands)
    if print_cmd:
        print(" ".join(commands))
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


def test_batomspy_location():
    batomspy = shutil.which("batomspy")
    assert batomspy is not None


def test_shebang():
    with tempfile.NamedTemporaryFile(suffix=".py", delete=False) as py_file:
        with open(py_file.name, "w") as fd:
            fd.write(BATOMSPY_TEST)
        os.chmod(py_file.name, 0o755)
    proc = _run_process([py_file.name])
    assert proc.returncode == 0


if __name__ == "__main__":
    test_batomspy_location()
    test_shebang()
    print("\n Shebang test: All pass! \n")
