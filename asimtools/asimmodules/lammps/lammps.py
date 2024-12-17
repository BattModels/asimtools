#!/usr/bin/env python
'''
Runs a user defined lammps script or template. LAMMPS must be installed 

Author: mkphuthi@github.com
'''
from typing import Dict, Optional
import sys
from pathlib import Path
from numpy.random import randint
import subprocess
import logging
from asimtools.utils import (
    get_atoms,
)

def lammps(
    template: str,
    image: Optional[Dict] = None,
    atom_style: str = 'atomic',
    variables: Optional[Dict] = None,
    placeholders: Optional[Dict] = None,
    lmp_cmd: str = 'lmp',
    masses: bool = True,
    seed: Optional[int] = None,
) -> Dict:
    """Runs a lammps script based on a specified template, variables can be 
    specified as arguments to be defined in the final LAMMPS input file if 
    placeholders are put in the template

    :param template: path to lammps input template file
    :type template: str
    :param image: Initial image for MD simulation. Image specification, 
        see :func:`asimtools.utils.get_atoms`, defaults to None
    :type image: Dict, optional
    :param atom_style: LAMMPS style in which to write image to Lammps data 
        input e.g. full, atomic etc., defaults to 'atomic'
    :type atom_style: str, optional
    :param variables: Dictionary of variables to be defined into the lammps 
        input, defaults to None
    :type variables: Dict, optional
    :param placeholders: Dictionary of placeholders to be filled into the 
        lammps input, defaults to None
    :type placeholders: Dict, optional
    :param lmp_cmd: Command with which to run LAMMPS, defaults to 'lmp'
    :type lmp_cmd: str, optional
    :param masses: Whether to specify atomic masses in LAMMPS data input, 
        requires ASE>3.23.0, defaults to True
    :type masses: bool, optional
    :param seed: Random seed for anywhere necessary in the template. You will 
        need to put the 'SEED'  placeholder anywhere you want a random 
        seed to be placed, if seed=None, a random one is generated, 
        defaults to None
    :type seed: int, optional
    :return: LAMMPS out file names
    :rtype: Dict
    """
    if variables is None:
        variables = {}

    # Make sure the provided asimmodule follows standard for reading
    # in arbitrary image provided by asimtools
    if image is not None:
        atoms = get_atoms(**image)
        if masses:
            try:
                atoms.write(
                    'image_input.lmpdat',
                    format='lammps-data',
                    atom_style=atom_style,
                    masses=masses,
                )
            except TypeError as te:
                err_txt = 'Need ASE version >=3.23 to support writing '
                err_txt += 'masses to lammps input file. Add mass keyword to '
                err_txt += 'lammps template instead'

                print(err_txt, file=sys.stderr)
                logging.error(err_txt)
                raise Exception(err_txt) from te

        else:
            atoms.write(
            'image_input.lmpdat',
            format='lammps-data',
            atom_style=atom_style,
        )
        variables['IMAGE_FILE'] = 'image_input.lmpdat'

    lmp_txt = ''
    for variable, value in variables.items():
        lmp_txt += f'variable {variable} equal {value}\n'

    lmp_txt += '\n'
    template = Path(template).resolve()
    with open(template, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    if placeholders is None:
        placeholders = {}

    for line in lines:
        if 'SEED' in line and 'SEED' not in placeholders:
            if seed is None:
                seed = str(randint(0, 100000))
            line = line.replace('SEED', seed)

        for placeholder in placeholders:
            line = line.replace(placeholder, placeholders[placeholder])
        lmp_txt += line

    if image is not None:
        assert 'read_data ' in lmp_txt, \
            'Make sure "read_data image_input.lmpdat" command is used \
                (with correct atom style) appropriately in lammps input \
                file if you specify image keyword'

    lmp_inp_file = 'input.lammps'
    with open(lmp_inp_file, 'w', encoding='utf-8') as f:
        f.write(lmp_txt)

    command = lmp_cmd + f' -i {lmp_inp_file}'
    command = command.split(' ')
    completed_process = subprocess.run(
            command, check=False, capture_output=True, text=True,
        )

    with open('lmp_stdout.txt', 'a+', encoding='utf-8') as f:
        f.write(completed_process.stdout)

    if completed_process.returncode != 0:
        err_txt = f'Failed to run {lmp_inp_file}\n'
        err_txt += 'See lmp_stderr.txt for details.'
        logging.error(err_txt)
        with open('lmp_stderr.txt', 'a+', encoding='utf-8') as f:
            f.write(completed_process.stderr)
        completed_process.check_returncode()
        return {}

    results = {'files': {
        'log': 'log.lammps',
        'thermo': 'log.lammps',
    }}

    return results
