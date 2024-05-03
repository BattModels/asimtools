from setuptools import setup, find_packages
from pathlib import Path
version_dict = {}
with open(Path(__file__).parents[0] / 'asimtools/_version.py') as fp:
    exec(fp.read(), version_dict)

setup(
    name="asimtools",
    version=version_dict["__version__"],
    description="Atomic Simulation Tools",
    author="Keith Phuthi and Emil Annevelink",
    python_requires=">=3.7",
    packages=find_packages(include=[
        "asimtools",
        "asimtools.*",
        "asimtools.scripts.*.*",
    ]),
    entry_points={
        # make the scripts available as command line scripts
        "console_scripts": [
            "asim-run = asimtools.scripts.asim_run:main",
            "asim-execute = asimtools.scripts.asim_execute:main",
            "asim-check = asimtools.scripts.asim_check:main",
        ]
    },
    install_requires=[
        "pandas",
        "pyyaml",
        "pymatgen",
        "ase>=3.22.1", 
        # Recommended to use the master branch of the ASE gitlab but only 
        # required for lammps to handle masses correctly as far as we can tell
        "colorama",
        "myst-parser",
        "sphinx",
        "pytest",
        "phonopy",
        "seekpath",
        "mp-api",
    ],
    zip_safe=True,
)
