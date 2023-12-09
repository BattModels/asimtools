#!/usr/bin/env python
'''
Apply the same asimmodule using multiple calculators

Author: mkphuthi@github.com

'''

from typing import Dict, Sequence, Optional, Union
from copy import deepcopy
from asimtools.job import DistributedJob
from asimtools.utils import get_calc_input, change_dict_value

def calc_array(
    subasimmodule_input: Dict,
    calc_ids: Sequence[str] = None,
    template_calc_id: Optional[str] = None,
    key_sequence: Optional[Sequence[str]] = None,
    array_values: Optional[Sequence] = None,
    env_ids: Optional[Union[Sequence[str],str]] = None,
    calc_input: Optional[Dict] = None,
    env_input: Optional[Dict] = None,
    ids: Sequence = None,
) -> Dict:
    """Apply the same asimmodule using different calculators and if necessary
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
    using_array_values = key_sequence is not None\
        and template_calc_id is not None\
        and array_values is not None
    assert calc_ids is not None or using_array_values, \
        'Specify either a sequence of calc_ids or all of key_sequence, \
        template_calc_id and array_values to iterate over'

    if using_array_values:
        if calc_input is None:
            calc_input = get_calc_input()
        calc_params = calc_input[template_calc_id]
        new_calc_input = {}
        if ids is None:
            ids = [f'{key_sequence[-1]}-{val}' for val in array_values]
        for i, val in enumerate(array_values):
            new_calc_params = change_dict_value(
                d=calc_params,
                new_value=val,
                key_sequence=key_sequence,
                return_copy=True,
            )
            new_calc_input[ids[i]] = new_calc_params
        calc_input = new_calc_input
        calc_ids = ids
    elif calc_ids is not None:
        if ids is None:
            ids = calc_ids

    assert len(ids) == len(calc_ids), \
        'Num. of calc_ids or array_values must match num. of ids'

    if env_ids is not None:
        assert len(env_ids) == len(calc_ids) or isinstance(env_ids, str), \
            'Provide one env_id or as many as there are calc_ids/array_values'

    array_sim_input = {}

    # Make individual sim_inputs for each calc
    for i, calc_id in enumerate(calc_ids):
        new_subasimmodule_input = deepcopy(subasimmodule_input)
        new_subasimmodule_input['args']['calc_id'] = calc_id
        array_sim_input[f'{ids[i]}'] = new_subasimmodule_input

        if env_ids is not None:
            if isinstance(env_ids, str):
                new_subasimmodule_input['env_id'] = env_ids
            else:
                new_subasimmodule_input['env_id'] = env_ids[i]

    # Create a distributed job object
    djob = DistributedJob(array_sim_input, env_input, calc_input)
    job_ids = djob.submit()

    results = {'job_ids': job_ids}
    return results
