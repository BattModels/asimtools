#!/usr/bin/env python
'''
Submits a chain of scripts in the order specified by steps
in input

Author: mkphuthi@github.com

'''

from typing import Dict, Tuple
from asimtools.job import branch, ChainedJob

@branch
def chained(
    config_input: Dict,
    steps: Dict,
    **kwargs
) -> Tuple[list,Dict]:
    '''
    Submits a chain of scripts in the order specified by steps
    in input
    '''

    cjob = ChainedJob(config_input, steps)
    job_ids = cjob.submit()

    results = cjob.unitjobs[-1].get_output()
    return job_ids, results
