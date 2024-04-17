from typing import Dict, Tuple, Sequence, Optional
from glob import glob
import os
from ase.io import read
from asimtools.calculators import load_calc
from asimtools.asimmodules.workflows.image_array import image_array
from asimtools.utils import get_str_btn, get_images

def forces(
    images: Dict,
    calc_id: str,
    calc_env_id: Optional[str] = None,
) -> Dict:

    calc = load_calc(calc_id)
    images = get_images(**images)
    singlepoint_input={
        'asimmodule': 'singlepoint',
        'args': {
            'calc_id': calc_id,
            'properties': ['energy', 'forces']
        },
    }
    if calc_env_id is not None:
        singlepoint_input.update({'env_id': calc_env_id,})

    results = image_array(
        images={
            'images': images,
        },
        labels=[f'supercell-{j:03d}' for j in range(len(images))],
        subsim_input=singlepoint_input,
    )

    return results
