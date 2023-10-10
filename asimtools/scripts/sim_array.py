#!/usr/bin/env python
'''
Runs the same script, iterating over multiple values of a specified argument
based on a sim_input template provided by the user

Author: mkphuthi@github.com

'''

from typing import Dict, Sequence, Optional, Union
from asimtools.job import DistributedJob
from asimtools.utils import change_dict_value

def sim_array(
    template_sim_input: Dict,
    key_sequence: Optional[Sequence[str]] = None,
    array_values: Optional[Sequence] = None,
    env_ids: Optional[Union[Sequence[str],str]] = None,
    calc_input: Optional[Dict] = None,
    env_input: Optional[Dict] = None,
    ids: Sequence = None,
) -> Dict:
    """Runs the same script, iterating over multiple values of a specified
    argument based on a sim_input template provided by the user

    :param template_sim_input: sim_input containing all the default parameters
    :type template_sim_input: Dict
    :param key_sequence: Sequence of keys to access value to be iterated over, defaults to None
    :type key_sequence: Optional[Sequence[str]], optional
    :param array_values: values to be iterated over in each simulation, defaults to None
    :type array_values: Optional[Sequence], optional
    :param env_ids: Environment(s) to be used for each simulation, must either be a list with as many env_ids as array values or a string with the env_id to be used by all simulations, defaults to None
    :type env_ids: Optional[Union[Sequence[str],str]], optional
    :param calc_input: calc_input file to use, defaults to None
    :type calc_input: Optional[Dict], optional
    :param env_input: env_input to use, defaults to None
    :type env_input: Optional[Dict], optional
    :param ids: Custom ids to use for each simulation, defaults to None
    :type ids: Sequence, optional
    :return: Results
    :rtype: Dict
    """
    if ids is None:
        ids = [f'{key_sequence[-1]}-{val}' for val in array_values]

    assert len(ids) == len(array_values), \
        'Num. of array_values must match num. of ids'

    if env_ids is None:
        pass
    elif isinstance(env_ids, str):
        env_ids = [env_ids for i in range(len(ids))]
    else:
        assert len(env_ids) == len(ids) or isinstance(env_ids, str), \
            'Provide one env_id or as many as there are array_values'

    # Create sim_inputs for each array value and to create input for
    # a DistributedJob object
    sim_inputs = {}
    for i, val in enumerate(array_values):
        new_sim_input = change_dict_value(
            d=template_sim_input,
            new_value=val,
            key_sequence=key_sequence,
            return_copy=True,
        )
        if env_ids is not None:
            new_sim_input['env_id'] = env_ids[i]
        sim_inputs[ids[i]] = new_sim_input

    # Create a distributed job object
    djob = DistributedJob(sim_inputs, env_input, calc_input)
    job_ids = djob.submit()

    results = {'job_ids': job_ids}
    return results
