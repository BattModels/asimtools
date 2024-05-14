'''
Global fixtures used in all tests defined here
'''

import pytest

@pytest.fixture
def inline_env_input():
    """Inline job env_input"""
    env_input = {
        'inline': {
            'mode': {
                'use_slurm': False,
                'interactive': True,
            },
        }
    }
    return env_input

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
def batch_dict_env_input():
    """Batch job env_input"""
    env_input = {
        'batch': {
            'mode': {
                'use_slurm': True,
                'interactive': False,
            },
            'slurm': {
                'flags': {'-n': 1, '-J': 'test', '--mem': '1G'}
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
def srun_env_input():
    """Interactive slurm job env_input"""
    env_input = {
        'srun': {
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
def all_env_input():
    """env_input with all test envs"""
    env_input = {
        'srun': {
            'mode': {
                'use_slurm': True,
                'interactive': True,
            },
            'slurm': {
                'flags': ['-n 1', '-J test', '--mem=1G']
            }
        },
        'batch': {
            'mode': {
                'use_slurm': True,
                'interactive': False,
            },
            'slurm': {
                'flags': {'-n': 1, '-J': 'test', '--mem': '1G'}
            }
        },
        'inline': {
            'mode': {
                'use_slurm': False,
                'interactive': True,
            },
        },
        'inline2': {
            'mode': {
                'use_slurm': False,
                'interactive': True,
            },
        },
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
        'asimmodule': 'singlepoint',
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
    Sim imput for a asimmodule that just completes without doing anything
    '''
    sim_input = {
        'asimmodule': './data/do_nothing.py',
        'prefix': 'test_',
        'args': {
            'duration': 3,
        },
    }
    return sim_input

@pytest.fixture
def lj_distributed_sim_input():
    ''' 
    Sim input for a distributed job that does some lj calculations
    '''
    subsim_input = {
        'asimmodule': 'singlepoint',
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
        'asimmodule': 'workflows.distributed',
        'env_id': 'inline',
        'args': {
            'subsim_inputs': {
                'id-0000': subsim_input,
                'id-0001': subsim_input,
                'id-0002': subsim_input,
            }
        }
    }

    return sim_input

@pytest.fixture
def lj_distributed_custom_name_sim_input():
    ''' 
    Sim input for a distributed job that does some lj calculations
    '''
    subsim_input = {
        'asimmodule': 'singlepoint',
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
        'asimmodule': 'workflows.distributed',
        'env_id': 'inline',
        'args': {
            'subsim_inputs': {
                'first': subsim_input,
                'second': subsim_input,
                'third': subsim_input,
                'id-03': subsim_input,
            }
        }
    }

    return sim_input

@pytest.fixture
def lj_distributed_batch_sim_input():
    ''' 
    Sim input for a distributed job that does some lj calculations
    '''
    subsim_input = {
        'asimmodule': 'singlepoint',
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
        'asimmodule': 'workflows.distributed',
        'env_id': 'batch',
        'args': {
            'subsim_inputs': {
                'id-0000': subsim_input,
                'id-0001': subsim_input,
                'id-0002': subsim_input,
                'id-0003': subsim_input,
            }
        }
    }

    return sim_input
