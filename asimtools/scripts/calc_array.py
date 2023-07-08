#!/usr/bin/env python
'''
Apply the same script using multiple calculators

Author: mkphuthi@github.com

'''

from typing import Dict, Sequence, Optional
from copy import deepcopy
from asimtools.job import DistributedJob

def calc_array(
    calc_ids: Sequence[str],
    subscript_input: Dict,
    env_ids: Optional[Sequence[str]] = None,
    calc_input: Optional[Dict] = None, #This doesn't work yet
    env_input: Optional[Dict] = None,
    ids: Sequence = None,
) -> Dict:
    """Apply the same script using different calculators and if necessary
    different environments

    :param calc_ids: Iterable with calc_ids, defaults to None
    :type calc_ids: Sequence, optional
    :param calc_input: Dictionary of calculator inputs
    :type calc_input: Dictionary, optional
    :param ids: Iterable with custom ids for each calc, defaults to None
    :type ids: Sequence, optional
    :return: Dictionary of results
    :rtype: Dict
    """
    array_sim_input = {}

    # Allow user to customize subdirectory names if needed
    if ids is None:
        ids = [calc_id for calc_id in calc_ids]
    else:
        assert len(ids) == len(calc_ids), \
            'Num. of calc_ids must match num. of ids'

    if env_ids is not None:
        assert len(env_ids) == len(calc_ids), \
            'Num. of calc_ids must match num. of env_ids'

    # Make individual sim_inputs for each calc
    for i, calc_id in enumerate(calc_ids):
        new_subscript_input = deepcopy(subscript_input)
        new_subscript_input['args']['calc_id'] = calc_id
        array_sim_input[f'{ids[i]}'] = new_subscript_input

        if env_ids is not None:
            new_subscript_input['env_id'] = env_ids[i]

    # Create a distributed job object
    djob = DistributedJob(array_sim_input, env_input, calc_input)
    job_ids = djob.submit()

    results = {'job_ids': job_ids}
    return results
