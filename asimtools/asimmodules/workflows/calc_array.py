#!/usr/bin/env python
'''
Apply the same asimmodule using multiple calculators

Author: mkphuthi@github.com

'''

from typing import Dict, Sequence, Optional, Union
from glob import glob
from copy import deepcopy
from asimtools.job import DistributedJob
from asimtools.utils import (
    get_calc_input,
    change_dict_value,
    get_str_btn,
)

def calc_array(
    subasimmodule_input: Dict,
    calc_ids: Sequence[str] = None,
    template_calc_id: Optional[str] = None,
    key_sequence: Optional[Sequence[str]] = None,
    array_values: Optional[Sequence] = None,
    env_ids: Optional[Union[Sequence[str],str]] = None,
    calc_input: Optional[Dict] = None,
    env_input: Optional[Dict] = None,
    file_pattern: Optional[str] = None,
    labels: Optional[Union[Sequence,str]] = 'values',
    str_btn_args: Optional[Dict] = None,
) -> Dict:
    """Apply the same asimmodule using different calculators and if necessary
    different environments

    :param calc_ids: Iterable with calc_ids, defaults to None
    :type calc_ids: Sequence, optional
    :param calc_input: Dictionary of calculator inputs
    :type calc_input: Dictionary, optional
    :param labels: Iterable with custom labels for each calc, defaults to None
    :type labels: Sequence, optional
    :return: Dictionary of results
    :rtype: Dict
    """
    using_array_values = key_sequence is not None\
        and template_calc_id is not None\
        and array_values is not None
    assert calc_ids is not None or using_array_values, \
        'Specify either an iterable of calc_ids or all of "key_sequence", \
        "template_calc_id" and "array_values" to iterate over'

    if using_array_values:
        if calc_input is None:
            calc_input = get_calc_input()
        calc_params = calc_input[template_calc_id]
        new_calc_input = {}
        if file_pattern is not None:
            array_values = glob(file_pattern)

        assert len(array_values) > 0, 'No array values or files found'
        if labels == 'str_btn':
            assert str_btn_args is not None, 'Provide str_btn_args for labels'
            labels = [get_str_btn(s, *str_btn_args) for s in array_values]
        elif labels == 'values':
            labels = [f'{key_sequence[-1]}-{val}' for val in array_values]
        elif labels is None:
            labels = [str(i) for i in range(len(array_values))]

        for i, val in enumerate(array_values):
            new_calc_params = change_dict_value(
                d=calc_params,
                new_value=val,
                key_sequence=key_sequence,
                return_copy=True,
            )
            new_calc_input[labels[i]] = new_calc_params
        calc_input = new_calc_input
        calc_ids = labels
    elif calc_ids is not None:
        assert labels != 'get_str_btn', \
            'get_str_btn only works when using the key_sequence argument.'
        if labels is None or labels == 'values':
            labels = calc_ids

    assert len(labels) == len(calc_ids), \
        'Num. of calc_ids or array_values must match num. of labels'

    if env_ids is not None:
        assert len(env_ids) == len(calc_ids) or isinstance(env_ids, str), \
            'Provide one env_id or as many as there are calc_ids/array_values'

    array_sim_input = {}

    # Make individual sim_inputs for each calc
    for i, calc_id in enumerate(calc_ids):
        new_subasimmodule_input = deepcopy(subasimmodule_input)
        new_subasimmodule_input['args']['calc_id'] = calc_id
        array_sim_input[f'{labels[i]}'] = new_subasimmodule_input

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
