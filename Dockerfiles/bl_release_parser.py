"""Parse the download link for a given blender module
"""
import os
import argparse
import requests
from bs4 import BeautifulSoup
from distutils.version import LooseVersion
import re
from subprocess import run
from pathlib import Path


DEFAULT_BLENDER_MIRROR = "https://mirrors.ocf.berkeley.edu/blender/release/"
BLENDER_MIRROR_ROOT = os.environ.get("BLENDER_MIRROR_URL", DEFAULT_BLENDER_MIRROR)

def get_latest_url(version="3.1", regex=r"blender-(\d+\.\d+\.\d+)-linux-x64.tar.xz"):
    url = f"{BLENDER_MIRROR_ROOT}/Blender{version}"
    req = requests.get(url)
    text = req.text
    soup = BeautifulSoup(text, "html.parser")
    hrefs = {}
    versions = []
    for link in soup.find_all("a"):
        href = link.get("href", "")
        if re.match(regex, href):
            match = next(re.finditer(regex, href))
            match_version = match.groups()[0]
            hrefs[match_version] = href
            versions.append(match_version)
    latest_href = hrefs[sorted(versions)[-1]]
    latest_url = f"{url}/{latest_href}"
    return latest_url, latest_href

def extract_blender(url, filename, root="/bin"):
    """Extract the tar.xz
    """
    commands = f"wget {url} && tar -xvf {filename} --strip-components=1 -C {root} && rm -rf ./blender-*"
    proc = run(commands, shell=True)
    if proc.returncode != 0:
        raise RuntimeError("Error extracting blender")
    

def get_default_bl_py_version(blender_root):
    """Get default python path
    """
    py_binary = next(Path(f"{blender_root}/python/bin/").glob("python*"))
    proc = run([py_binary.as_posix(), "-V"], capture_output=True)
    version = proc.stdout.decode("ascii").strip().split()[-1]
    return version

def extract_python_source(blender_root, python_version):
    short_version = ".".join(python_version.split(".")[:2])
    url = f"https://www.python.org/ftp/python/{python_version}/Python-{python_version}.tgz"
    commands = f"wget url && tar -xzf Python-*.tgz && cp -r Python-*/Include/* {blender_root}/python/include/python{short_version}/ && rm -rf Python-*"
    proc = run(commands, shell=True)
    if proc.returncode != 0:
        raise RuntimeError("Error extracting python source")



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("version", help="Blender version in major.minor split")
    arguments = parser.parse_args()
    version = arguments.version
    url, href = get_latest_url(version)
    print(url, href)
    extract_blender(url, href, root="/bin")
    blender_root = f"/bin/{version}"
    python_version = get_default_bl_py_version()
    extract_python_source(blender_root, python_version)
    return



if __name__ == "__main__":
    main()

