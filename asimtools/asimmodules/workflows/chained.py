#!/usr/bin/env python
'''
Submits a chain of asimmodules in the order specified by steps
in input

Author: mkphuthi@github.com

'''

from typing import Dict, Optional
from asimtools.job import ChainedJob


def chained(
    steps: Dict,
    env_input: Optional[Dict] = None,
    calc_input: Optional[Dict] = None,
) -> Dict:
    """Run a workflow using the provided sim_input.yaml one after the other,
    The jobs will run one after the other, if slurm is used, job dependencies
    will be used. Note that if one of the subsim_inputs internally launches another
    asimmodule e.g. image_array, calc_array etc, slurm job dependencies will not
    work as expected. In those cases, we recommend running the steps one by one
    or interactively


    :param steps: Dictionary where each key is an ID and the value is a \
        standard sim_input
    :type steps: Dict
    :param env_input: env_input that overrides global file, defaults to None
    :type env_input: Optional[Dict], optional
    :param calc_input: calc_input that overrides global file, defaults to None
    :type calc_input: Optional[Dict], optional
    :return: Dictionary of results
    :rtype: Dict
    """
    cjob = ChainedJob(
        sim_input=steps, env_input=env_input, calc_input=calc_input
    )
    job_ids = cjob.submit()

    results = cjob.unitjobs[-1].get_output()
    results.update({'job_ids': job_ids})
    return results
