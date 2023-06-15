#!/usr/bin/env python
'''
Submits a chain of scripts in the order specified by steps
in input

Author: mkphuthi@github.com

'''

from typing import Dict, Tuple, Union
from asimtools.job import ChainedJob


def chained(
    steps: Dict,
    env_input: Union[Dict,None] = None,
    calc_input: Union[Dict,None] = None,
) -> Tuple[list,Dict]:
    '''
    Submits a chain of scripts in the order specified by steps
    in input
    '''
    cjob = ChainedJob(
        sim_input=steps, env_input=env_input, calc_input=calc_input
    )
    job_ids = cjob.submit()

    results = cjob.unitjobs[-1].get_output()
    results.update({'job_ids': job_ids})
    return results
