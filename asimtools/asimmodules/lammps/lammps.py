#!/usr/bin/env python
'''
Runs a user defined lammps script or template 

Author: mkphuthi@github.com
'''
from typing import Dict
from pathlib import Path
import subprocess
import logging
from asimtools.utils import (
    get_atoms,
    join_names,
)

def lammps(
    template: str,
    image: Dict = None,
    prefix: str = '',
    atom_style: str = 'atomic',
    variables: Dict = None,
    lmp_cmd: str = 'lmp',
    masses: bool = True,
) -> Dict:
    ''' 
    Runs a lammps simulation based on a template lammps input script
    '''
    if variables is None:
        variables = {}

    # Make sure the provided asimmodule follows standard for reading
    # in arbitrary image provided by asimtools
    if image is not None:
        atoms = get_atoms(**image)
        if masses:
            atoms.write(
                'image_input.lmpdat',
                format='lammps-data',
                atom_style=atom_style,
                masses=masses,
            )
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

    for line in lines:
        lmp_txt += line

    if image is not None:
        assert 'read_data ' in lmp_txt, \
            'Make sure "read_data image_input.lmpdat" command is used \
                (with correct atom style) appropriately in lammps input \
                asimmodule if you specify image keyword'

    lmp_inp_file = join_names([prefix, 'input.lammps'])
    with open(lmp_inp_file, 'w', encoding='utf-8') as f:
        f.write(lmp_txt)

    command = lmp_cmd + f' -i {lmp_inp_file}'
    command = command.split(' ')
    completed_process = subprocess.run(
            command, check=False, capture_output=True, text=True,
        )

    with open('lmp_stdout.txt', 'w', encoding='utf-8') as f:
        f.write(completed_process.stdout)

    if completed_process.returncode != 0:
        err_txt = f'Failed to run {lmp_inp_file}\n'
        err_txt += 'See lmp.stderr.txt for details.'
        logging.error(err_txt)
        with open('lmp_stderr.txt', 'w', encoding='utf-8') as f:
            f.write(completed_process.stderr)
        completed_process.check_returncode()
        return {}

    results = {'files': {
        'log': 'log.lammps',
        'thermo': 'log.lammps',
    }}

    return results