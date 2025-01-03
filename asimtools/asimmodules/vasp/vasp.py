#!/usr/bin/env python
'''
Runs VASP based on input files and optionally MP settings. 
Heavily uses pymatgen for IO and MP settings.
VASP must be installed 

Author: mkphuthi@github.com
'''
from typing import Dict, Optional, Sequence
import os
import sys
from pathlib import Path
from numpy.random import randint
import subprocess
import logging
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp.sets import (
    MPRelaxSet, MPStaticSet, MPNonSCFSet, MPScanRelaxSet
)
from asimtools.utils import (
    get_atoms,
)

def vasp(
    image: Optional[Dict],
    vaspinput_args: Optional[Dict] = None,
    command: str = 'vasp_std',
    mpset: Optional[str] = None,
    write_image_output: bool = True,
) -> Dict:
    """Run VASP with given input files and specified image

    :param image: Initial image for VASP calculation. Image specification,
        see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param vaspinput_args: Dictionary of pymatgen's VASPInput arguments. 
        See :class:`pymatgen.io.vasp.inputs.VaspInput`
    :type vaspinput_args: Dict
    :param command: Command with which to run VASP, defaults to 'vasp_std'
    :type command: str, optional
    :param mpset: Materials Project VASP set to use, defaults to None
    :type mpset: str, optional
    :param write_image_output: Whether to write output image in standard 
        asimtools format to file, defaults to False
    :type write_image_output: bool, optional
    """

    if vaspinput_args:
        if image is not None:
            atoms = get_atoms(**image)
            struct = AseAtomsAdapter.get_structure(atoms)
        vaspinput = VaspInput(
            poscar=Poscar(struct),
            **vaspinput_args
        )
    else:
        atoms = get_atoms(**image)
        struct = AseAtomsAdaptor.get_structure(atoms)
        if mpset == 'MPRelaxSet':
            vasp_input = MPRelaxSet(struct)
        elif mpset == 'MPStaticSet':
            vasp_input = MPStaticSet(struct)
        elif mpset == 'MPNonSCFSet':
            vasp_input = MPNonSCFSet(struct)
        elif mpset == 'MPScanRelaxSet':
            vasp_input = MPScanRelaxSet(struct)
        else:
            raise ValueError(f'Unknown MPSet: {mpset}')

        vasp_input.write_input("./")

    command = command.split(' ')
    completed_process = subprocess.run(
            command, check=False, capture_output=True, text=True,
        )

    with open('vasp_stdout.txt', 'a+', encoding='utf-8') as f:
        f.write(completed_process.stdout)

    if completed_process.returncode != 0:
        err_txt = f'Failed to run VASP\n'
        err_txt += 'See vasp_stderr.txt for details.'
        logging.error(err_txt)
        with open('vasp_stderr.txt', 'a+', encoding='utf-8') as f:
            f.write(completed_process.stderr)
        completed_process.check_returncode()
        return {}

    if write_image_output:
        atoms_output = read('OUTCAR')
        atoms_output.write(
            'image_output.xyz',
            format='extxyz',
        )

    return {}
