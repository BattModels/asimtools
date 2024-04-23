#!/usr/bin/env python
'''
Calculates single point energy 

Author: mkphuthi@github.com
'''

from typing import Tuple, Dict, Optional
import logging
from asimtools.calculators import load_calc
from asimtools.utils import (
    get_atoms,
)

def singlepoint(
    calc_id: str,
    image: Dict,
    properties: Tuple[str] = ('energy', 'forces'),
    prefix: Optional[str] = None,
) -> Dict:
    """Evaluates the properties of a single image, currently implemented
    properties are energy, forces and stress

    :param calc_id: calc_id specification
    :type calc_id: str
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param properties: properties to evaluate, defaults to ('energy', 'forces')
    :type properties: Tuple[str], optional
    :return: Dictionary of results
    :rtype: Dict
    """
    calc = load_calc(calc_id)
    atoms = get_atoms(**image)
    atoms.set_calculator(calc)

    if prefix is not None:
        prefix = prefix + '_'
    else:
        prefix = ''
    import time; time.sleep(20)
    if 'energy' in properties:
        try:
            energy = atoms.get_potential_energy()
            logging.debug('Calculated energy')
        except Exception:
            logging.error('Failed to calculate energy')
            raise

    if 'forces' in properties:
        try:
            atoms.get_forces()
            logging.debug('Calculated forces')
        except Exception:
            logging.error('Failed to calculate forces')
            raise

    if 'stress' in properties:
        try:
            atoms.get_stress()
            logging.debug('Calculated stress')
        except Exception:
            logging.error('Failed to calculate stress')
            raise

    image_file = prefix + 'image_output.xyz'
    atoms.write(image_file, format='extxyz')

    results = {
        'energy': float(energy),
        'files': {'image': image_file}
    }
    return results
