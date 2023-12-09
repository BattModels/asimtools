#!/usr/bin/env python
'''
Describe the asimmodule briefly here. If this is a asimmodule that runs multiple steps,
describe it here using reStructuredText to generate autodocs

Cite the papers where the method/asimmodule was first introduced here as well

Author: mkphuthi@github.com
'''

from typing import Dict, Sequence, Optional
import logging
from asimtools.calculators import load_calc
from asimtools.utils import (
    get_atoms,
)

def vacancy_formation_energy(
    calc_id: str,
    image: Dict,
) -> Dict:
    '''

    '''
    calc = load_calc(calc_id)
    atoms = get_atoms(**image)
    atoms.set_calculator(calc)

    results = {}
    return results
