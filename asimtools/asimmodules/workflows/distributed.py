#!/usr/bin/env python
'''
Distributes a number of simulations at the same time if possible

Author: mkphuthi@github.com

'''

from typing import Dict, Optional
from asimtools.job import DistributedJob

def distributed(
    subsim_inputs: Dict,
    env_input: Optional[Dict] = None,
    calc_input: Optional[Dict] = None,
    array_max: Optional[int] = None,
) -> Dict:
    """Distributes a set of jobs based on inpout parameter. The scnarios are
    as follows:
    #. use_slurm: False, interactive: True -> Runs one after the other.
    #. use_slurm: True, interactive: False, same env_id -> Launches a job \
    # array.
    #. use_slurm: False, interactive: False, different env_ids -> Launches \
    # multiple jobs.

    :param subsim_inputs: Dictionary of asimmodules, each key is an ID and each \
        value is a sim_input file
    :type subsim_inputs: Dict
    :param env_input: env_input that overrides global file, defaults to None
    :type env_input: Optional[Dict], optional
    :param calc_input: calc_input that overrides global file, defaults to None
    :type calc_input: Optional[Dict], optional
    :param array_max: Maximum number of jobs to run in array, defaults to None
    :type array_max: Optional[int], optional
    :return: Dictionary of results
    :rtype: Dict
    """
    djob = DistributedJob(
        sim_input=subsim_inputs,
        env_input=env_input,
        calc_input=calc_input
    )
    job_ids = djob.submit(array_max=array_max)

    results = {'job_ids': job_ids}
    return results
