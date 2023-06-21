#!/usr/bin/env python
'''
Distributes a number of simulations at the same time if possible

Author: mkphuthi@github.com

'''

from typing import Dict, Optional
from asimtools.job import DistributedJob

def distributed(
    subscripts: Dict,
    env_input: Optional[Dict] = None,
    calc_input: Optional[Dict] = None,
) -> Dict:
    '''
    Submits multiple scripts simultaneously and handles
    file management based on standard
    '''

    djob = DistributedJob(
        sim_input=subscripts,
        env_input=env_input,
        calc_input=calc_input
    )
    job_ids = djob.submit()

    results = {'job_ids': job_ids}
    return results
