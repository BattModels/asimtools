#!/usr/bin/env python
'''
Distributes a number of simulations at the same level

Author: mkphuthi@github.com

'''

from typing import Dict
from asimtools.job import leaf, ChainedJob

@leaf
def chained(
    config_input: Dict,
    steps: Dict,
    **kwargs
):
    '''
    Submits multiple scripts simultaneously and handles
    file management based on standard
    '''

    cjob = ChainedJob(config_input, steps)
    cjob.submit()

    results = cjob.unitjobs[-1].get_output()
    return results
