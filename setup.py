from setuptools import setup, find_packages

setup(
    name='batoms',
    version='2.3.0',
    packages=find_packages(),
    license='GPL',
    description='A Python package for creating, editing and rendering atoms and molecules structures using Blender.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Beautiful Atoms Team',
    author_email='xing.wang@gmail.com',
    url='https://github.com/beautiful-atoms/beautiful-atoms',
    install_requires=[
        "bpy",
        "ase",
    ],
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        'Programming Language :: Python :: 3',
    ],
    python_requires='>=3.10',
    entry_points={
    },
)
