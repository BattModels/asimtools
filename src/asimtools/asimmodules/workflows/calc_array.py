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
    calculators: Sequence[Dict] = None,
    template_calculator: Optional[Dict] = None,
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
    skip_failed: Optional[bool] = False,
    group_size: int = 1,
) -> Dict:
    """Apply the same asimmodule using different calculators and if necessary
    different environments

    :param calc_ids: Deprecated. Use calculators instead. Iterable with
        calc_ids, defaults to None
    :type calc_ids: Sequence, optional
    :param template_calc_id: Deprecated. Use template_calculator instead.
        calc_id of the template calculator, defaults to None
    :type template_calc_id: str, optional
    :param calculators: Sequence of calculator dicts, each with a 'calc_id'
        or 'calc_params' key, see :func:`asimtools.calculators.load_calc`,
        defaults to None
    :type calculators: Sequence[Dict], optional
    :param template_calculator: Calculator dict with 'calc_id' or 'calc_params'
        key used as template when iterating over key_sequence values, see
        :func:`asimtools.calculators.load_calc`, defaults to None
    :type template_calculator: Optional[Dict], optional
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
    :param array_max: Maximum number of jobs to run simultanteously in job
        array, defaults to None
    :type array_max: Optional[int], optional
    :param skip_failed: Skip failed jobs and move to next, defaults to False
    :type skip_failed: Optional[bool], optional
    :param group_size: Number of jobs to group together, defaults to 1
    :type group_size: int, optional
    :return: Dictionary of results
    :rtype: Dict
    """
    # Backward compatibility: convert old-style params to new interface
    if calculators is None and calc_ids is not None:
        calculators = [{'calc_id': cid} for cid in calc_ids]
    if template_calculator is None and template_calc_id is not None:
        template_calculator = {'calc_id': template_calc_id}

    print([
            array_values, linspace_args, arange_args, file_pattern
        ])
    using_array_values = key_sequence is not None\
        and template_calculator is not None\
        and [
            array_values, linspace_args, arange_args, file_pattern
        ].count(None) == 3
    err_txt = 'Specify either a sequence of "calculators" or all of '
    err_txt += '"key_sequence", "template_calculator" and one of ['
    err_txt += '"array_values", "linspace_args", "arange_args", "file_pattern"'
    err_txt += '] to iterate over'
    assert calculators is not None or using_array_values, err_txt

    if using_array_values:
        if template_calculator.get('calc_params') is not None:
            calc_params = template_calculator['calc_params']
        else:
            if calc_input is None:
                calc_input = get_calc_input()
            calc_params = calc_input[template_calculator['calc_id']]

        if calculators is not None and labels is None:
            labels = [c.get('calc_id', f'calc-{i}') for i, c in enumerate(calculators)]

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

        calculators = []
        for i, val in enumerate(array_values):
            new_calc_params = change_dict_value(
                dct=calc_params,
                new_value=val,
                key_sequence=key_sequence,
                return_copy=True,
            )
            if secondary_array_values is not None:
                for k, vs in zip(secondary_key_sequences, secondary_array_values):
                    new_calc_params = change_dict_value(
                        dct=new_calc_params,
                        new_value=vs[i],
                        key_sequence=k,
                        return_copy=False,
                    )
            calculators.append({'calc_params': new_calc_params})

    elif calculators is not None:
        assert labels != 'get_str_btn', \
            'get_str_btn only works when using the key_sequence argument.'
        if labels is None or labels == 'values':
            labels = [c.get('calc_id', f'calc-{i}') for i, c in enumerate(calculators)]

    assert len(labels) == len(calculators), \
        'Num. of calculators or array_values must match num. of labels'

    if env_ids is not None:
        assert len(env_ids) == len(calculators) or isinstance(env_ids, str), \
            'Provide one env_id or as many as there are calculators/array_values'

    array_sim_input = {}

    # Make individual sim_inputs for each calc
    for i, calculator in enumerate(calculators):
        new_subsim_input = deepcopy(subsim_input)
        new_subsim_input['args']['calculator'] = calculator
        array_sim_input[f'{labels[i]}'] = new_subsim_input

        if env_ids is not None:
            if isinstance(env_ids, str):
                new_subsim_input['env_id'] = env_ids
            else:
                new_subsim_input['env_id'] = env_ids[i]

    # Create a distributed job object
    djob = DistributedJob(array_sim_input, env_input, calc_input)
    job_ids = djob.submit(
        array_max=array_max,
        skip_failed=skip_failed,
        group_size=group_size,
    )

    results = {'job_ids': job_ids}
    return results
