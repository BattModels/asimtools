#!/usr/bin/env python
'''
Apply the same script on multiple images

Author: mkphuthi@github.com

'''

from typing import Dict, Sequence
from copy import deepcopy
from asimtools.job import leaf, DistributedJob
from asimtools.utils import get_images

@leaf
def image_array(
    config_input: Dict,
    images: Dict,
    subscript_input: Dict,
    ids: Sequence[str] = None,
    **kwargs
):
    ''' Submits same script on multiple images '''
    images = get_images(**images)
    array_sim_input = {}

    # Allow user to customize subdirectories if needed
    if ids is None:
        ids = [str(i) for i in range(len(images))]
    else:
        assert len(ids) == len(images), \
            'Num. of images must match nim. of ids'

    for i, image in enumerate(images):
        new_subscript_input = deepcopy(subscript_input)
        new_subscript_input['args']['image'] = {
            'atoms': image
        }
        array_sim_input[f'{ids[i]}'] = new_subscript_input

    djob = DistributedJob(config_input, array_sim_input)
    djob.submit()

    results = {}
    return results
