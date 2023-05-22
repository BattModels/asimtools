'''.Xauthority'''

from typing import Dict, Tuple
import numpy as np
from asimtools.job import branch, UnitJob
from asimtools.utils import (
    get_atoms,
)

@branch
def singlepoint_calculations(
    config_input: Dict,
    image: Dict,
    singlepoint_config_id: str,
    preprocess_config_id: str,
    nimages: int = 5,
    scale_range: Tuple[float,float] = (0.95, 1.05),
    **kwargs
):
    ''' Does calculations for eos '''

    assert scale_range[0] < scale_range[1], 'x_scale[0] should be smaller than\
         x_scale[1]'

    atoms = get_atoms(**image)

    # Preprocessing the data, structures and data structures
    scales = np.linspace(scale_range[0], scale_range[1], nimages)
    scaled_atoms = []
    scaled_ids = []
    for scale in scales:
        new_cell = atoms.get_cell() * scale
        new_atoms = atoms.copy()
        new_atoms.set_cell(new_cell)
        scaled_atoms.append(new_atoms)
        scaled_ids.append(f'x{scale:.2f}')

    image_array_input = {
        'script': 'image_array',
        'config_id': preprocess_config_id,
        'args': {
            'workdir': '.',
            'images': scaled_atoms,
            'ids': scaled_ids,
            'subscript_input': {
                'script': 'singlepoint',
                'config_id': singlepoint_config_id,
                'args': {
                    'properties': ['energy'],
                }
            }
        }
    }

    image_array_job = UnitJob(config_input, image_array_input)
    image_array_job.submit()
    return {}
