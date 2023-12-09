#!/usr/bin/env python
'''
Describe the asimmodule briefly here. If this is a asimmodule that runs multiple steps,
describe it here using reStructuredText to generate autodocs

Cite the papers where the method/asimmodule was first introduced here as well

Author: mkphuthi@github.com
'''

from typing import Dict 
# from asimtools.calculators import load_calc
from asimtools.utils import (
    get_atoms,
    # get_images,
)

def template(
    # calc_id: str,
    # image: Dict, # Optional: if using a structure
    # images: Dict, # Optional: if using multiple structures
) -> Dict:
    '''
    asimmodule does xyz specifically
    '''

    # calc = load_calc(calc_id) # Optional: Any code that loads a calculator is replaced by this one line!
    # atoms = get_atoms(**image) # Optional
    # images = get_images(**images) # Optional
    # atoms.set_calculator(calc) # Optional

    #############################
    # Do some cool science here #
    #############################

    results = {
        'energy': 'Placeholder_for_my_result',
        'files': {'image': 'placeholder_four_image_output'}
    }
    return results # Always return a dictionary! Use {} if necessary
