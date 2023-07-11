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
    """Run the same script using multiple env_ids

    :param subscript_input: sim_input for the script to be run
    :type subscript_input: Dict
    :param env_ids: Sequence of env_ids to iterate over
    :type env_ids: Sequence[str]
    :param calc_input: calc_input to override global file, defaults to None
    :type calc_input: Optional[Dict], optional
    :param env_input: env_input defaults to overide global file, defaults to\
        None
    :type env_input: Optional[Dict], optional
    :param ids: Sequence of custom ids for each env_id, defaults to None
    :type ids: Sequence, optional
    :return: Dictionary of results
    :rtype: Dict
    """
    array_sim_input = {}

    # Allow user to customize subdirectory names if needed
    if ids is None:
        ids = [str(i) for i in range(len(env_ids))]
    else:
        assert len(ids) == len(env_ids), \
            'Num. of env_ids must match num. of ids'

    # Make individual sim_inputs for each env
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
