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
    **kwargs,
) -> Dict:
    """Calculate forces from displaced images for phonopy calculations

    :param images: Images specification, see :func:`asimtools.utils.get_images`
    :type images: Dict
    :param calc_id: calc_id specification, see :func:`asimtools.utils.get_calc`
    :type calc_id: str
    :param calc_env_id: env_id to use for the calculator, defaults to None
    :type calc_env_id: Optional[str], optional
    :param kwargs: Additional keyword arguments to pass to image_array
    :type kwargs: Dict
    :return: Results
    :rtype: Dict
    """
    images = get_images(**images)
    singlepoint_input={
        'asimmodule': 'singlepoint',
        'args': {
            'calc_id': calc_id,
            'properties': ['energy', 'forces']
        },
    }

    results = image_array(
        images={
            'images': images,
        },
        labels=[f'supercell-{j:03d}' for j in range(len(images))],
        subsim_input=singlepoint_input,
        env_ids=calc_env_id,
        **kwargs
    )

    return results
