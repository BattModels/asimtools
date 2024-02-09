#!/usr/bin/env python
'''
Runs the same asimmodule, iterating over multiple values of a specified argument
based on a sim_input template provided by the user

Author: mkphuthi@github.com

'''

from typing import Dict, Sequence, Optional, Union
from glob import glob
from copy import deepcopy
import numpy as np
from asimtools.job import DistributedJob
from asimtools.utils import change_dict_value
from asimtools.utils import get_str_btn

def sim_array(
    template_sim_input: Dict,
    key_sequence: Optional[Sequence[str]] = None,
    array_values: Optional[Sequence] = None,
    file_pattern: Optional[str] = None,
    linspace_args: Optional[Sequence] = None,
    arange_args: Optional[Sequence] = None,
    env_ids: Optional[Union[Sequence[str],str]] = None,
    calc_input: Optional[Dict] = None,
    env_input: Optional[Dict] = None,
    labels: Optional[Union[Sequence,str]] = 'values',
    str_btn_args: Optional[Dict] = None,
    secondary_key_sequences: Optional[Sequence] = None,
    secondary_array_values: Optional[Sequence] = None,
) -> Dict:
    """Runs the same asimmodule, iterating over multiple values of a specified
    argument based on a sim_input template provided by the user

    :param template_sim_input: sim_input containing all the default parameters
    :type template_sim_input: Dict
    :param key_sequence: Sequence of keys to access value to be iterated over,
        defaults to None
    :type key_sequence: Optional[Sequence[str]], optional
    :param array_values: values to be iterated over in each simulation,
        defaults to None
    :type array_values: Optional[Sequence], optional
    :param file_pattern: pattern of files to be iterated over in each
        simulation, defaults to None
    :type st: Optional[str], optional
    :param linspace_args: arguments to pass to :func:`numpy.linspace` to be
        iterated over in each simulation, defaults to None
    :type linspace_args: Optional[Sequence], optional
    :param arange_args: arguments to pass to :func:`numpy.arange` to be
        iterated over in each simulation, defaults to None
    :type arange_args: Optional[Sequence], optional
    :param labels: Custom labels to use for each simulation, defaults to None
    :type labels: Sequence, optional
    :param env_ids: Environment(s) to be used for each simulation, must either
        be a list with as many env_ids as array values or a string with the
        env_id to be used by all simulations, defaults to None
    :type env_ids: Optional[Union[Sequence[str],str]], optional
    :param secondary_key_sequences: list of other keys to iterate over in
        tandem with key_sequence to allow changing multiple key-value pairs,
        defaults to None
    :type secondary_key_sequences: Sequence, optional
    :param secondary_array_values: list of other other array_values to iterate
        over in tandem with array_values to allow changing multiple key-value
        pairs, defaults to None
    :type secondary_array_values: Sequence, optional
    :param calc_input: calc_input file to use, defaults to None
    :type calc_input: Optional[Dict], optional
    :param env_input: env_input to use, defaults to None
    :type env_input: Optional[Dict], optional
    :return: Results
    :rtype: Dict
    """

    possible_values = [array_values, file_pattern, linspace_args, arange_args]
    assert possible_values.count(None) + 1 == len(possible_values), \
        f'Provide only one of {possible_values}'

    if file_pattern is not None:
        array_values = glob(file_pattern)
    elif linspace_args is not None:
        array_values = np.linspace(*linspace_args)
    elif arange_args is not None:
        array_values = np.arange(*arange_args)

    assert len(array_values) > 0, 'No array values or files found'

    if labels == 'str_btn':
        assert str_btn_args is not None, 'Provide str_btn_args for labels'
        labels = [get_str_btn(s, *str_btn_args) for s in array_values]
    elif labels == 'values':
        if key_sequence is not None:
            labels = [f'{key_sequence[-1]}-{val}' for val in array_values]
        else:
            labels = [f'value-{val}' for val in array_values]
    elif labels is None:
        labels = [str(i) for i in range(len(array_values))]

    assert len(labels) == len(array_values), \
        'Num. of array_values must match num. of labels'

    if secondary_array_values is not None:
        for l in secondary_array_values:
            assert len(l) == len(labels), \
                f"Secondary values ({len(l)}) not same length as array values"\
                f" ({len(labels)})"

    if env_ids is None:
        pass
    elif isinstance(env_ids, str):
        env_ids = [env_ids for i in range(len(labels))]
    else:
        assert len(env_ids) == len(labels) or isinstance(env_ids, str), \
            'Provide one env_id or as many as there are array_values'

    # Create sim_inputs for each array value and to create input for
    # a DistributedJob object
    sim_inputs = {}
    for i, val in enumerate(array_values):
        if key_sequence is not None:
            new_sim_input = change_dict_value(
                d=template_sim_input,
                new_value=val,
                key_sequence=key_sequence,
                return_copy=True,
            )
        else:
            new_sim_input = deepcopy(template_sim_input)

        if secondary_array_values is not None:
            for k, vs in zip(secondary_key_sequences, secondary_array_values):
                print('kv', k, vs[i])
                new_sim_input = change_dict_value(
                    d=new_sim_input,
                    new_value=vs[i],
                    key_sequence=k,
                    return_copy=False,
                )

        if env_ids is not None:
            new_sim_input['env_id'] = env_ids[i]
        sim_inputs[labels[i]] = new_sim_input

    # Create a distributed job object
    djob = DistributedJob(sim_inputs, env_input, calc_input)
    job_ids = djob.submit()

    results = {'job_ids': job_ids}
    return results
