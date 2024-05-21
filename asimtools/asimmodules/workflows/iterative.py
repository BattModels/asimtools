#!/usr/bin/env python
'''
Runs the same asimmodule, iterating over multiple values of a specified
argument one after the other based on a sim_input template provided by the user

Author: mkphuthi@github.com

'''

from typing import Dict, Sequence, Optional, Union
import os
from pathlib import Path
from copy import deepcopy
from asimtools.job import ChainedJob
from asimtools.utils import change_dict_value
from asimtools.asimmodules.workflows.utils import prepare_array_vals

def iterative(
    template_sim_input: Dict,
    dependent_file: Optional[os.PathLike] = None,
    dependent_file_key_sequence: Optional[Sequence[str]] = None,
    key_sequence: Optional[Sequence[str]] = None,
    array_values: Optional[Sequence] = None,
    file_pattern: Optional[str] = None,
    linspace_args: Optional[Sequence] = None,
    arange_args: Optional[Sequence] = None,
    env_ids: Optional[Union[Sequence[str],str]] = None,
    calc_input: Optional[Dict] = None,
    env_input: Optional[Dict] = None,
    # labels: Optional[Union[Sequence,str]] = 'values',
    # label_prefix: Optional[str] = None,
    str_btn_args: Optional[Dict] = None,
    secondary_key_sequences: Optional[Sequence] = None,
    secondary_array_values: Optional[Sequence] = None,
    # array_max: Optional[int] = None,
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
    :type file_pattern: Optional[str], optional
    :param linspace_args: arguments to pass to :func:`numpy.linspace` to be
        iterated over in each simulation, defaults to None
    :type linspace_args: Optional[Sequence], optional
    :param arange_args: arguments to pass to :func:`numpy.arange` to be
        iterated over in each simulation, defaults to None
    :type arange_args: Optional[Sequence], optional
    :param labels: Custom labels to use for each simulation, defaults to None
    :type labels: Sequence, optional
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
    :param array_max: Number of jobs to run at once in scheduler
    :type array_max: int, optional
    :param calc_input: calc_input file to use, defaults to None
    :type calc_input: Optional[Dict], optional
    :param env_input: env_input to use, defaults to None
    :type env_input: Optional[Dict], optional
    :return: Results
    :rtype: Dict
    """

    results = prepare_array_vals(
        key_sequence=key_sequence,
        array_values=array_values,
        file_pattern=file_pattern,
        linspace_args=linspace_args,
        arange_args=arange_args,
        env_ids=env_ids,
        # labels=None,
        # label_prefix=label_prefix,
        str_btn_args=str_btn_args,
        secondary_key_sequences=secondary_key_sequences,
        secondary_array_values=secondary_array_values,
    )
    array_values = results['array_values']
    print(array_values)
    labels = [f'step-{i}' for i in range(len(array_values))]
    env_ids = results['env_ids']
    secondary_array_values = results['secondary_array_values']
    secondary_key_sequences = results['secondary_key_sequences']

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

        if dependent_file_key_sequence is not None and i > 0:
            dep_arg = str(Path(f'../step-{i-1}') / dependent_file)
            new_sim_input = change_dict_value(
                d=new_sim_input,
                new_value=dep_arg,
                key_sequence=dependent_file_key_sequence,
                return_copy=False,
            )
        else:
            new_sim_input = deepcopy(new_sim_input)

        if secondary_array_values is not None:
            for k, vs in zip(secondary_key_sequences, secondary_array_values):
                new_sim_input = change_dict_value(
                    d=new_sim_input,
                    new_value=vs[i],
                    key_sequence=k,
                    return_copy=False,
                )

        if env_ids is not None:
            new_sim_input['env_id'] = env_ids[i]
        sim_inputs[labels[i]] = new_sim_input

    # Create a chained job object
    cjob = ChainedJob(sim_inputs, env_input, calc_input)
    job_ids = cjob.submit()

    results = {'job_ids': job_ids}
    return results
