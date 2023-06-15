#!/usr/bin/env python
'''
Distributes a number of simulations at the same time if possible

Author: mkphuthi@github.com

'''

from typing import Dict, Tuple
from asimtools.job import branch, DistributedJob

@branch
def distributed(
    config_input: Dict,
    subscripts: Dict,
    **kwargs
) -> Tuple[list,Dict]:
    '''
    Submits multiple scripts simultaneously and handles
    file management based on standard
    '''

    djob = DistributedJob(config_input, subscripts)
    job_ids = djob.submit()

    results = {'job_ids': job_ids}
    return results
