#!/usr/bin/env python
'''
Apply the same script using different envs

Author: mkphuthi@github.com

'''

from typing import Dict, Sequence, Optional
from copy import deepcopy
from asimtools.job import DistributedJob

def env_array(
    subscript_input: Dict,
    env_ids: Sequence[str],
    calc_input: Optional[Dict] = None,
    env_input: Optional[Dict] = None,
    ids: Sequence = None,
) -> Dict:
    ''' Submits same script on multiple images '''
    array_sim_input = {}

    # Allow user to customize subdirectory names if needed
    if ids is None:
        ids = [str(i) for i in range(len(env_ids))]
    else:
        assert len(ids) == len(env_ids), \
            'Num. of images must match num. of ids'

    # Make individual sim_inputs for each image
    for i, env_id in enumerate(env_ids):
        new_subscript_input = deepcopy(subscript_input)
        new_subscript_input['env_id'] = env_id
        array_sim_input[f'{ids[i]}'] = new_subscript_input

    # Create a distributed job object
    djob = DistributedJob(array_sim_input, env_input, calc_input)
    job_ids = djob.submit()

    results = {
        'job_ids': job_ids,
    }
    return results
