'''
Global fixtures used in all tests defined here
'''

import pytest

@pytest.fixture
def batch_env_input():
    """Batch job env_input"""
    env_input = {
        'batch': {
            'mode': {
                'use_slurm': True,
                'interactive': False,
            },
            'slurm': {
                'flags': ['-n 1', '-J test', '--mem=1G']
            }
        }
    }
    return env_input

@pytest.fixture
def salloc_env_input():
    """Interactive slurm job env_input"""
    env_input = {
        'salloc': {
            'mode': {
                'use_slurm': True,
                'interactive': True,
            },
            'slurm': {
                'flags': ['-n 1', '-J test', '--mem=1G']
            }
        }
    }
    return env_input

@pytest.fixture
def inline_env_input():
    """Inline job env_input"""
    env_input = {
        'batch': {
            'mode': {
                'use_slurm': False,
                'interactive': True,
            },
        }
    }
    return env_input

@pytest.fixture
def lj_argon_calc_input():
    """ Calc input that uses LennardJones in a batch job"""
    calc_input = {
        'lj': {
            'name': 'LennardJones',
            'module': 'ase.calculators.lj',
            'args': {
                'sigma': 3.54,
                'epsilon': 0.00802236,
            },
        },
    }
    return calc_input

@pytest.fixture
def singlepoint_argon_sim_input():
    ''' Sim input for singlepoint calculation '''
    sim_input = {
        'script': 'singlepoint',
        'prefix': 'test_',
        'image': {
            'name': 'Ar',
            'crystalstructure': 'fcc',
            'a': 5.4,
        },
        'properties': ['energy', 'force', 'stress'],
    }
    return sim_input

@pytest.fixture
def do_nothing_sim_input():
    ''' 
    Sim imput for a script that just completes without doing anything
    '''
    sim_input = {
        'script': './data/do_nothing.py',
        'prefix': 'test_',
        'args': {
            'duration': 5,
        },
    }
    return sim_input
