'''.Xauthority'''

from typing import Dict, Tuple
import numpy as np
from ase.io import write
from asimtools.utils import (
    get_atoms,
)

def preprocess(
    image: Dict,
    nimages: int = 5,
    scale_range: Tuple[float,float] = (0.95, 1.05),
) -> Dict:
    ''' Prepares structures for EOS calculation '''

    assert scale_range[0] < scale_range[1], 'x_scale[0] should be smaller than\
         x_scale[1]'

    atoms = get_atoms(**image)

    # Make a database of structures with the volumes scaled appropriately
    scales = np.linspace(scale_range[0], scale_range[1], nimages)
    scaled_images = []
    for scale in scales:
        new_cell = atoms.get_cell() * scale
        new_atoms = atoms.copy()
        new_atoms.set_cell(new_cell)
        new_atoms.info['scale'] = f'{scale:.2f}'
        scaled_images.append(new_atoms)

    scaled_images_file = 'preprocess_images_output.xyz'
    write(scaled_images_file, scaled_images, format='extxyz')

    return {'files':{'images': scaled_images_file}}
