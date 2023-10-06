#!/usr/bin/env python
'''
Runs a user defined lammps script or template 

Author: mkphuthi@github.com
'''
from typing import Dict, Sequence
from asimtools.job import UnitJob
from copy import deepcopy

def entropy_from_covariance(
    image: Dict,
    temps: Sequence,
    calc_env_id: str,
    repeat: Sequence[int],
    process_env_id: str,
    md_script: str,
    md_args: Dict,
    sampling_index: str = ':',
    system: str = 'mod',
) -> Dict:
    ''' 
    Runs the entropy calculation using xdat package on a given image for the 
    given temperatures
    '''

    md_args['image'] = deepcopy(image)
    md_args['image']['repeat'] = repeat
    sim_input = {
        'script': 'distributed',
        'env_id': process_env_id,
        'submit': True,
        'overwrite': True,
        'args': {
            'subscripts': {},
        }
    }
    print('image', image)
    for temp in temps:
        if md_script == 'lammps':
            md_args['variables']['TEMP'] = temp
        else:
            md_args['temp'] = temp
        temp_input = {
            'script': 'chained',
            'env_id': process_env_id,
            'submit': True,
            'overwrite': True,
            'args': {
                'steps': {
                    'step-0': {
                        'script': md_script,
                        'env_id': calc_env_id,
                        'submit': True,
                        'overwrite': True,
                        'args': deepcopy(md_args),
                    },
                    'step-1': {
                        'script': 'entropy_from_covariance.run_xdat',
                        'env_id': process_env_id,
                        'submit': True,
                        'overwrite': True,
                        'args': {
                            'unit_cell_image': image,
                            'ideal_image': {
                                'image_file': '../step-0/input_image.xyz',
                            },
                            'md_images': {
                                'image_file': '../step-0/output.traj',
                                'index': sampling_index,
                            },
                            'temperature': temp,
                            'system': system,
                        },
                    },
                },
            },
        }

        sim_input['args']['subscripts'][f'T{temp:.1f}K'] = temp_input
        sim_input['overwrite'] = True
    results = UnitJob(sim_input).submit()

    return results