from setuptools import setup, find_packages
# from pathlib import Path
# see https://packaging.python.org/guides/single-sourcing-package-version/
# version_dict = {}
# with open(Path(__file__).parents[0] / "nequip/_version.py") as fp:
#     exec(fp.read(), version_dict)
# version = version_dict["__version__"]
# del version_dict

setup(
    name="asimtools",
    description="Atomic Simulation Tools",
    author="Keith Phuthi and Emil Annevelink",
    python_requires=">=3.9",
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
        "ase",
    ],
    zip_safe=True,
)
