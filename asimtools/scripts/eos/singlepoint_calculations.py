'''.Xauthority'''

from typing import Dict, Tuple
import numpy as np
# from asimtools.job import branch
from asimtools.scripts.image_array import image_array
from asimtools.utils import (
    get_atoms,
)

# @branch
def singlepoint_calculations(
    singlepoint_env_id: str,
    calc_input: str,
    calc_id: str,
    env_id: str,
    image: Dict,
    nimages: int = 5,
    scale_range: Tuple[float,float] = (0.95, 1.05),
    **kwargs
) -> Tuple[None,Dict]:
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

    subscript_input = {
        'script': 'singlepoint',
        'env_id': singlepoint_env_id,
        'calc_params': calc_input[calc_id],
        'args': {
            'properties': ['energy'],
        }
    }
    results = image_array(
        images={'images': scaled_atoms},
        subscript_input=subscript_input,
        env_input=None,
        calc_input=None,
        ids=scaled_ids
    )
    return results
