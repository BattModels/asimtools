#!/usr/bin/env python
'''
Calculates single point energy 

Author: mkphuthi@github.com
'''

#pylint: disable=unused-argument
# pylint: disable=too-many-arguments

from typing import Tuple, Dict
from asimtools.calculators import load_calc
from asimtools.job import Job
from asimtools.utils import (
    get_atoms,
    join_names,
)

def singlepoint(
    calc_input: Dict,
    image: Dict = None,
    prefix: str = '',
    properties: Tuple[str] = ('energy', 'forces'),
    **kwargs
) -> Dict:
    ''' 
    Calculates the single point energy, forces and stresses where possible
    '''
    sim_input = {'prefix': prefix}
    job = Job(calc_input, sim_input)
    job.start()

    calc = load_calc(calc_input)

    atoms = get_atoms(**image)
    atoms.set_calculator(calc)

    if 'energy' in properties:
        try:
            energy = atoms.get_potential_energy()
        except Exception:
            job.update_status('failed')
            print('Failed to calculate energy')
            raise

    if 'forces' in properties:
        try:
            atoms.get_forces()
        except Exception:
            job.update_status('failed')
            print('Failed to calculate forces')
            raise

    if 'stress' in properties:
        try:
            atoms.get_stress()
        except Exception:
            job.update_status('failed')
            print('Failed to calculate stress')
            raise

    image_file = join_names([prefix, 'image_output.xyz'])
    atoms.write(image_file, format='extxyz')

    job.update_output({'energy': float(energy)})
    job.add_output_files({'image': image_file})

    job.complete()
    return job.get_output()
