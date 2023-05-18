'''
Global fixtures used in all tests defined here
'''

import pytest

@pytest.fixture
def slurm_batch_calc_input():
    ''' Calc input that uses no calculator in a batch job '''
    calc_input = {
        'job': {
            'use_slurm': True,
            'interactive': False,
        },
        'slurm': {
            'flags': ['-n 2', '--gres=gpu:1']
        }
    }
    return calc_input

@pytest.fixture
def lj_slurm_batch_calc_input():
    ''' Calc input that uses LennardJones in a batch job '''
    calc_input = {
        'name': 'LennardJones',
        'job': {
            'use_slurm': True,
            'interactive': False,
        },
        'slurm': {
            'flags': ['-n 2', '--gres=gpu:1']
        }
    }
    return calc_input

@pytest.fixture
def interactive_calc_input():
    ''' Calc input that uses no calculator in the terminal '''
    calc_input = {
        'job': {
            'use_slurm': False,
            'interactive': True,
        }
    }
    return calc_input

@pytest.fixture
def lj_interactive_calc_input():
    ''' Calc input that uses LennardJones in the terminal '''
    calc_input = {
        'name': 'LennardJones',
        'job': {
            'use_slurm': False,
            'interactive': True,
        }
    }
    return calc_input

@pytest.fixture
def singlepoint_sim_input():
    ''' Sim input for singlepoint calculation '''
    sim_input = {
        'script': 'singlepoint',
        'prefix': 'test_',
        'image': {
            'name': 'Ar',
            'crystalstructure': 'fcc',
            'a': 5.4,
        }
    }
    return sim_input

@pytest.fixture
def do_nothing_sim_input():
    ''' 
    Sim imput for a script that just completes without doing anything
    '''
    sim_input = {
        'script': './data/do_nothing_script.py',
        'prefix': 'test_',
    }
    return sim_input
