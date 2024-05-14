#!/usr/bin/env python
'''
Describe the asimmodule briefly here. If this is a asimmodule that runs multiple steps,
describe it here using reStructuredText to generate autodocs

Cite the papers where the method/asimmodule was first introduced here as well

Author: mkphuthi@github.com
'''

from typing import Dict 
# from asimtools.calculators import load_calc
# from asimtools.utils import (
#     get_atoms,
#     get_images,
# )

def template(
    # calc_id: str,
    # image: Dict, # Optional: if using a structure
    # images: Dict, # Optional: if using multiple structures
) -> Dict:
    '''
    asimmodule does xyz
    '''

    # Some useful lines
    # calc = load_calc(calc_id)
    # atoms = get_atoms(**image)
    # images = get_images(**images)

    #############################
    # Do some cool science here #
    #############################

    results = {
        'my_result': 'placeholder_for_my_result',
        'files': {'file_description': 'placeholder_for_output_file'}
    }
    return results # Always return a dictionary! Use {} if necessary
