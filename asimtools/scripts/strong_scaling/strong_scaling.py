'''
Very broken, idea is to standarize how we do strong scaling for 
DFT calculation benchmarks
'''

from typing import Dict, Sequence
from copy import deepcopy
from asimtools.job import ChainedJob

# @branch
def strong_scaling(
    subscript_input: Dict,
    subscript_env_input: Dict,
    ntasks: Sequence[int],
    nnodes: Sequence[int] = None,
    ntask_flag: str = '-n',
    log_x: bool = True,
    log_y: bool = True,
) -> Dict:
    ''' Create env_array to measure strong scaling of a script '''
    assert subscript_env_input['mode'].get('use_slurm', False), \
        'Scaling only makes sense if "use_slurm=True" in subscript_env_input'
    assert subscript_env_input.get('slurm', {}).get('flags', False), \
        'Provide some slurm flags'

    flags = subscript_env_input.get('slurm', {}).get('flags', [])
    for flag in flags:
        assert '-n ' not in flag, f'Do not specify {flag} in slurm \
            parameters in subscript_env_input, use ntasks arg'
        assert '-N ' not in flag, f'Do not specify {flag} in slurm \
            parameters in subscript_env_input, use nnodes arg'

    if nnodes is None:
        nnodes = [1] * len(ntasks)
    else:
        ntasks_per_node = False
        for flag in flags:
            if '--ntasks-per-node' in flag:
                ntasks_per_node = int(flag.split('=')[-1])
        assert ntasks_per_node, \
            'If you scale nodes, specify --ntasks-per-node in env_input'

    ids = [f'n{nprocess}' for nprocess in ntasks]
    env_input = {}
    for i, nid in enumerate(ids):
        env_input[nid] = deepcopy(subscript_env_input)
        env_input[nid]['slurm']['flags'].append(
            f'{ntask_flag} {ntasks[i]}'
        )

        env_input[nid]['slurm']['flags'].append(
            f'-N {nnodes[i]}'
        )

    sim_input = {
        'step-0': {
            'script': 'env_array',
            'args': {
                'subscript_input': subscript_input,
                'ids': ids,
                'env_ids': ids,
                'env_input': env_input
            },
        },
        'step-1': {
            'script': 'strong_scaling.postprocess',
            'submit': False,
            'args': {
                'outdir_pattern': '../step-0/id-*',         
                'log_x': log_x,
                'log_y': log_y,
            },
        },
    }

    cjob = ChainedJob(sim_input, env_input=None, calc_input=None)
    cjob.submit()
    return cjob.get_output()
