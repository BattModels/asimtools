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

@pytest.fixture
def do_nothing_distributed_sim_input():
    ''' 
    Sim imput for a script that just completes without doing anything
    '''
    subscript_input = {
        'script': 'singlepoint',
        'env_id': 'inline',
        'args': {
            'calc_id': 'lj',
            'image': {
                'name': 'Ar',
            },
            'properties': ['energy', 'forces'],
        }
    }
    sim_input = {
        'script': 'distributed',
        'env_id': 'inline',
        'args': {
            'subscripts': {
                'id-00': subscript_input,
                'id-01': subscript_input,
                'id-02': subscript_input,
                'id-03': subscript_input,
                'id-04': subscript_input,
                'id-05': subscript_input,
                'id-06': subscript_input,
                'id-07': subscript_input,
                'id-08': subscript_input,
                'id-09': subscript_input,
                'id-10': subscript_input,
                'id-11': subscript_input,
                'id-12': subscript_input,
            }
        }
    }

    return sim_input

@pytest.fixture
def do_nothing_distributed_custom_name_sim_input():
    ''' 
    Sim imput for a script that just completes without doing anything
    '''
    subscript_input = {
        'script': 'singlepoint',
        'env_id': 'inline',
        'args': {
            'calc_id': 'lj',
            'image': {
                'name': 'Ar',
            },
            'properties': ['energy', 'forces'],
        }
    }
    sim_input = {
        'script': 'distributed',
        'env_id': 'inline',
        'args': {
            'subscripts': {
                'first': subscript_input,
                'second': subscript_input,
                'third': subscript_input,'id-03': subscript_input,
                'id-04': subscript_input,
                'id-05': subscript_input,
                'id-06': subscript_input,
                'id-07': subscript_input,
                'id-08': subscript_input,
                'id-09': subscript_input,
                'id-10': subscript_input,
                'id-11': subscript_input,
                'id-12': subscript_input,
            }
        }
    }

    return sim_input

@pytest.fixture
def do_nothing_distributed_batch_sim_input():
    ''' 
    Sim imput for a script that just completes without doing anything
    '''
    subscript_input = {
        'script': 'singlepoint',
        'env_id': 'inline',
        'args': {
            'calc_id': 'lj',
            'image': {
                'name': 'Ar',
            },
            'properties': ['energy', 'forces'],
        }
    }
    sim_input = {
        'script': 'distributed',
        'env_id': 'batch',
        'args': {
            'subscripts': {
                'id-00': subscript_input,
                'id-01': subscript_input,
                'id-02': subscript_input,
                'id-03': subscript_input,
                'id-04': subscript_input,
                'id-05': subscript_input,
                'id-06': subscript_input,
                'id-07': subscript_input,
                'id-08': subscript_input,
                'id-09': subscript_input,
                'id-10': subscript_input,
                'id-11': subscript_input,
                'id-12': subscript_input,
            }
        }
    }

    return sim_input
