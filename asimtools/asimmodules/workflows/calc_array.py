#!/usr/bin/env python
'''
Apply the same asimmodule using multiple calculators

Author: mkphuthi@github.com

'''

from typing import Dict, Sequence, Optional, Union
from copy import deepcopy
from asimtools.job import DistributedJob
from asimtools.utils import (
    get_calc_input,
    change_dict_value,
)
from asimtools.asimmodules.workflows.utils import prepare_array_vals

def calc_array(
    subsim_input: Dict,
    calc_ids: Sequence[str] = None,
    template_calc_id: Optional[str] = None,
    key_sequence: Optional[Sequence[str]] = None,
    array_values: Optional[Sequence] = None,
    file_pattern: Optional[str] = None,
    linspace_args: Optional[Sequence] = None,
    arange_args: Optional[Sequence] = None,
    label_prefix: Optional[str] = None,
    env_ids: Optional[Union[Sequence[str],str]] = None,
    calc_input: Optional[Dict] = None,
    env_input: Optional[Dict] = None,
    labels: Optional[Union[Sequence,str]] = 'values',
    str_btn_args: Optional[Dict] = None,
    secondary_key_sequences: Optional[Sequence] = None,
    secondary_array_values: Optional[Sequence] = None,
    array_max: Optional[int] = None,
) -> Dict:
    """Apply the same asimmodule using different calculators and if necessary
    different environments

    :param calc_ids: Iterable with calc_ids, defaults to None
    :type calc_ids: Sequence, optional
    :param calc_input: Dictionary of calculator inputs
    :type calc_input: Dictionary, optional
    :param labels: Iterable with custom labels for each calc, defaults to None
    :type labels: Sequence, optional
        :param array_values: values to be iterated over in each simulation,
        defaults to None
    :type array_values: Optional[Sequence], optional
    :param file_pattern: pattern of files to be iterated over in each
        simulation, defaults to None
    :type file_pattern: Optional[str], optional
    :param linspace_args: arguments to pass to :func:`numpy.linspace` to be
        iterated over in each simulation, defaults to None
    :type linspace_args: Optional[Sequence], optional
    :param arange_args: arguments to pass to :func:`numpy.arange` to be
        iterated over in each simulation, defaults to None
    :type arange_args: Optional[Sequence], optional
    :param label_prefix: Prefix to add before labels which can make extracting 
        data from file paths easier, defaults to None
    :type label_prefix: str, optional
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
    :return: Dictionary of results
    :rtype: Dict
    """
    print([
            array_values, linspace_args, arange_args, file_pattern
        ])
    using_array_values = key_sequence is not None\
        and template_calc_id is not None\
        and [
            array_values, linspace_args, arange_args, file_pattern
        ].count(None) == 3
    err_txt = 'Specify either a sequence of "calc_ids" or all of '
    err_txt += '"key_sequence", "template_calc_id" and one of ['
    err_txt += '"array_values", "linspace_args", "arange_args", "file_pattern"'
    err_txt += '] to iterate over'
    assert calc_ids is not None or using_array_values, err_txt

    if using_array_values:
        if calc_input is None:
            calc_input = get_calc_input()
        calc_params = calc_input[template_calc_id]
        new_calc_input = {}

        if calc_ids is not None and labels is None:
            labels = calc_ids

        results = prepare_array_vals(
            key_sequence=key_sequence,
            array_values=array_values,
            file_pattern=file_pattern,
            linspace_args=linspace_args,
            arange_args=arange_args,
            env_ids=env_ids,
            labels=labels,
            label_prefix=label_prefix,
            str_btn_args=str_btn_args,
            secondary_key_sequences=secondary_key_sequences,
            secondary_array_values=secondary_array_values,
        )
        array_values = results['array_values']
        labels = results['labels']
        env_ids = results['env_ids']
        secondary_array_values = results['secondary_array_values']
        secondary_key_sequences = results['secondary_key_sequences']

        for i, val in enumerate(array_values):
            new_calc_params = change_dict_value(
                d=calc_params,
                new_value=val,
                key_sequence=key_sequence,
                return_copy=True,
            )
            if secondary_array_values is not None:
                for k, vs in zip(secondary_key_sequences, secondary_array_values):
                    new_calc_params = change_dict_value(
                        d=new_calc_params,
                        new_value=vs[i],
                        key_sequence=k,
                        return_copy=False,
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
        new_subsim_input = deepcopy(subsim_input)
        new_subsim_input['args']['calc_id'] = calc_id
        array_sim_input[f'{labels[i]}'] = new_subsim_input

        if env_ids is not None:
            if isinstance(env_ids, str):
                new_subsim_input['env_id'] = env_ids
            else:
                new_subsim_input['env_id'] = env_ids[i]

    # Create a distributed job object
    djob = DistributedJob(array_sim_input, env_input, calc_input)
    job_ids = djob.submit(array_max=array_max)

    results = {'job_ids': job_ids}
    return results
