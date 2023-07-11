#!/usr/bin/env python
'''
Calculates single point energy 

Author: mkphuthi@github.com
'''

from typing import Tuple, Dict
from asimtools.calculators import load_calc
# from asimtools.job import uses_calc
from asimtools.utils import (
    get_atoms,
)

def singlepoint(
    calc_id: str,
    image: Dict,
    properties: Tuple[str] = ('energy', 'forces'),
) -> Dict:
    """Evaluates the properties of a single image, currently implemented
    properties are energy, forces and stress

    :param calc_id: ID of calculator provided in calc_input or global file
    :type calc_id: str
    :param image: Image config, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param properties: properties to evaluate, defaults to ('energy', 'forces')
    :type properties: Tuple[str], optional
    :return: Dictionary of results
    :rtype: Dict
    """
    calc = load_calc(calc_id)
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

    image_file = 'image_output.xyz'
    atoms.write(image_file, format='extxyz')

    results = {
        'energy': float(energy),
        'files': {'image': image_file}
    }
    return results
