#!/usr/bin/env python
'''
Runs a user defined lammps script or template 

### BROKEN ###

Author: mkphuthi@github.com
'''
from typing import Dict
import subprocess
from asimtools.utils import (
    get_atoms,
    join_names,
)

# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=dangerous-default-value

def lammps(
    template: str,
    image: Dict = None,
    prefix: str = '',
    atom_style: str = 'full',
    variables: Dict = {},
    lmp_cmd: str = 'lmp',
) -> Dict:
    ''' 
    Runs a lammps simulation based on a template lammps input script
    '''

    lmp_txt = ''
    for variable, value in variables.items():
        lmp_txt += f'variable {variable} equal {value}\n'

    lmp_txt += '\n'
    with open(template, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    for line in lines:
        lmp_txt += line

    # Make sure the provided script follows standard for reading
    # in arbitrary image provided by asimtools
    if image is not None:
        assert 'read_data ${image_file}' in lmp_txt, \
            'Make sure "read_data ${image_file}" command is used (with correct atom style) \
            in lammps input script if you specify image keyword'
        atoms = get_atoms(**image)
        atoms.write(
            'image_input.lmpdat', format='lammps-data', atom_style=atom_style
        )

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
        print(err_txt)
        with open('lmp_stderr.txt', 'w', encoding='utf-8') as f:
            f.write(completed_process.stderr)
        completed_process.check_returncode()
        return None

    results = {'files': {
        'log': 'log.lammps',
        'thermo': 'log.lammps',
        # 'dump': 'dump'
    }}

    return results
