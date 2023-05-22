#!/usr/bin/env python
'''
Distributes a number of simulations at the same level

Author: mkphuthi@github.com

'''

from typing import Dict
from asimtools.job import branch, DistributedJob

@branch
def distributed(
    config_input: Dict,
    subscripts: Dict,
    **kwargs
):
    '''
    Submits multiple scripts simultaneously and handles
    file management based on standard
    '''

    djob = DistributedJob(config_input, subscripts)
    djob.submit()

    results = {}
    return results
