#!/usr/bin/env python
'''
Calculates single point energy 

Author: mkphuthi@github.com
'''

# pylint: disable=unused-argument
# pylint: disable=too-many-arguments

from typing import Tuple, Dict
from asimtools.calculators import load_calc
from asimtools.job import leaf
from asimtools.utils import (
    get_atoms,
    join_names,
)

@leaf
def singlepoint(
    config_input: Dict,
    image: Dict = None,
    prefix: str = '',
    properties: Tuple[str] = ('energy', 'forces'),
    **kwargs
) -> Dict:
    ''' 
    Calculates the single point energy, forces and stresses where possible
    '''
    calc = load_calc(config_input['calc'])
    atoms = get_atoms(**image)

    atoms.set_calculator(calc)

    if 'energy' in properties:
        try:
            energy = atoms.get_potential_energy()
        except Exception:
            print('Failed to calculate energy')
            raise

    if 'forces' in properties:
        try:
            atoms.get_forces()
        except Exception:
            print('Failed to calculate forces')
            raise

    if 'stress' in properties:
        try:
            atoms.get_stress()
        except Exception:
            print('Failed to calculate stress')
            raise

    image_file = join_names([prefix, 'image_output.xyz'])
    atoms.write(image_file, format='extxyz')

    results = {
        'energy': float(energy),
        'files': {'image': image_file}
    }
    return results
