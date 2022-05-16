import sys
import subprocess
import importlib

dependencies = {"ase": "ase",
                "scikit-image": "skimage",
                }

def check_module(modname):
    try:
        module = importlib.import_module(modname)
        return True
    except ImportError:
        return False


def install_pip():
    # check pip exist or not
    # if not exist
    if subprocess.call([sys.executable, "-m", "pip", "--version"]):
        cmd = [sys.executable, "-m", "ensurepip", "--upgrade"]
        subprocess.call(cmd)
        # upgrade
        cmd = [sys.executable, "-m", "pip", "install", "--upgrade", "pip"]
        subprocess.call(cmd)
    else:
        print("pip exist.")


def install():
    install_pip()
    for package, modname in dependencies.items():
        if not check_module(modname):
            cmd = [sys.executable, "-m", "pip",
                   "install", "--upgrade", package]
            subprocess.call(cmd)
        else:
            print("package {0} installed.".format(package))
