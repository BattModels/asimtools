#!/usr/bin/env python
'''
Apply the same script on multiple images

Author: mkphuthi@github.com

'''

from typing import Dict, Optional, List
from copy import deepcopy
from asimtools.job import DistributedJob

def strong_scaling(
    ntask_array: List[int],
    subscript_sim_input: Dict,
    subscript_env_input: Optional[Dict] = None,
    calc_input: Optional[Dict] = None,
) -> Dict:
    ''' Submits same script on multiple images '''

    assert subscript_env_input['mode'].get('use_slurm', False), \
        'Scaling only makes sense if "use_slurm=True in subscript_env_input'

    flags = subscript_env_input.get('slurm', {}).get('flags', [])
    for flag in flags:
        assert '-n' not in flag, 'Do not specify tasks in slurm parameters in \
            subscript_env_input e.g. ' + {flag}

    env_inputs = {}
    array_sim_input = {}
    for ntask in ntask_array:
        env_input = deepcopy(subscript_env_input)
        if not env_input.get('slurm', False):
            env_input['slurm'] = {}
        if len(flags) == 0:
            env_input['slurm']['flags'] = []
        if not env_input['slurm'].get('precommands', False):
            env_input['slurm']['precommands'] = []
        if not env_input['slurm'].get('postcommands', False):
            env_input['slurm']['postcommands'] = []

        env_input['slurm']['flags'].append(f'-n {ntask}')
        env_input['slurm']['precommands'].append('start=$SECONDS')
        env_input['slurm']['postcommands'].append('duration=$(( SECONDS - start))')
        env_input['slurm']['postcommands'].append('echo')
        ntask_id = f'ntask-{ntask}'
        env_inputs[ntask_id] = env_input
        subscript_sim_input['env_id'] = ntask_id
        array_sim_input[ntask_id] = subscript_sim_input

    # Create a distributed job object
    djob = DistributedJob(array_sim_input, env_inputs, calc_input)
    job_ids = djob.submit()

    results = {'job_ids': job_ids}
    return results
