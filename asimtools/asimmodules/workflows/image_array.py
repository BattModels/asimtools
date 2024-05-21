#!/usr/bin/env python
'''
Apply the same asimmodule on multiple images

Author: mkphuthi@github.com

'''

from typing import Dict, Sequence, Optional, Union
from copy import deepcopy
from asimtools.job import DistributedJob
from asimtools.utils import get_images, change_dict_value
from asimtools.asimmodules.workflows.utils import prepare_array_vals

def image_array(
    images: Dict,
    subsim_input: Dict,
    calc_input: Optional[Dict] = None,
    env_input: Optional[Dict] = None,
    array_max: Optional[int] = None,
    key_sequence: Optional[Sequence[str]] = None,
    env_ids: Optional[Union[Sequence[str],str]] = None,
    labels: Optional[Union[Sequence,str]] = None,
    label_prefix: Optional[str] = None,
    str_btn_args: Optional[Dict] = None,
    secondary_key_sequences: Optional[Sequence] = None,
    secondary_array_values: Optional[Sequence] = None,
) -> Dict:
    """Submit the same asimmodule on multiple images and if specified, use
    different env_ids

    :param images: Images specification, see :func:`asimtools.utils.get_images`
    :type images: Dict
    :param subsim_input: sim_input of asimmodule to be run
    :type subsim_input: Dict
    :param calc_input: calc_input to override global file, defaults to None
    :type calc_input: Optional[Dict], optional
    :param env_input: env_input to override global file, defaults to None
    :type env_input: Optional[Dict], optional
    :param array_max: Maximum jobs to run simultanteously in job array, \
        defaults to None
    :type array_max: Optional[int], optional
    :param labels: Custom labels for each image, defaults to None
    :type labels: Sequence[str], optional
    :param env_ids: Sequence of envs for each image, must be the same length \
        as the total number of images, defaults to None
    :type env_ids: Sequence[str], optional
    :param labels: Custom labels to use for each simulation, defaults to \
        None
    :param label_prefix: Prefix to add before labels which can make extracting 
        data from file paths easier, defaults to None
    :type label_prefix: str, optional
    :param secondary_key_sequences: list of other keys to iterate over in
        tandem with key_sequence to allow changing multiple key-value pairs,
        defaults to None
    :type secondary_key_sequences: Sequence, optional
    :param secondary_array_values: list of other other array_values to iterate
        over in tandem with images to allow changing multiple key-value
        pairs, defaults to None
    :type secondary_array_values: Sequence, optional
    :return: Dictionary of results
    :rtype: Dict
    """
    images = get_images(**images)
    array_values = images

    if labels is None:
        labels = [str(i) for i in range(len(array_values))]
    results = prepare_array_vals(
        array_values=array_values,
        key_sequence=key_sequence,
        env_ids=env_ids,
        labels=labels,
        label_prefix=label_prefix,
        str_btn_args=str_btn_args,
        secondary_key_sequences=secondary_key_sequences,
        secondary_array_values=secondary_array_values,
    )

    if key_sequence is None:
        key_sequence = ['args', 'image']
        # For backwards compatibility where we don't have to specify image
        subsim_input = deepcopy(subsim_input)
        subsim_input['args']['image'] = {}
    key_sequence += ['atoms']

    labels = results['labels']
    env_ids = results['env_ids']
    secondary_array_values = results['secondary_array_values']
    secondary_key_sequences = results['secondary_key_sequences']

    # Create sim_inputs for each array value and to create input for
    # a DistributedJob object
    array_sim_input = {}
    for i, val in enumerate(array_values):
        new_sim_input = change_dict_value(
            d=subsim_input,
            new_value=val,
            key_sequence=key_sequence,
            return_copy=True,
        )

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

        array_sim_input[f'{labels[i]}'] = new_sim_input

    # Create a distributed job object
    djob = DistributedJob(array_sim_input, env_input, calc_input)
    job_ids = djob.submit(array_max=array_max)

    results = {'job_ids': job_ids}
    return results
