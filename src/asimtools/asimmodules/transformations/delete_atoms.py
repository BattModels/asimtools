'''
Produce a set of images with unit cells scaled compared to the input

author: mkphuthi@github.com
'''

from typing import Dict, Optional, Sequence
import numpy as np
from ase.io import write
from asimtools.utils import (
    get_atoms,
)

def apply_scale(old_atoms, scale):
    ''' Applies a scaling factor to a unit cell '''
    atoms = old_atoms.copy()
    new_cell = atoms.get_cell() * scale
    atoms.set_cell(new_cell, scale_atoms=True)
    atoms.info['scale'] = f'{scale:.3f}'
    return atoms

def scale_unit_cells(
    image: Dict,
    scales:  Optional[Sequence] = None,
    logspace: Optional[Sequence] = None,
    linspace: Optional[Sequence] = None,
    scale_by: str = 'a',
) -> Dict:
    """Produce a set of images with unit cells scaled compared to the input

    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param scales: Scaling values by which to scale cell, defaults to None
    :type scales: Optional[Sequence], optional
    :param logspace: Parameters to pass to np.logspace for scaling values,
        defaults to None
    :type logspace: Optional[Sequence], optional
    :param linspace: Parameters to pass to np.linspace for scaling values,
        defaults to None
    :type linspace: Optional[Sequence], optional
    :param scale_by: Scale either "volume" or "a" which is lattice parameter,
        defaults to 'a'
    :type scale_by: str, optional
    :raises ValueError: If more than one of scales, linspace, logspace are
        provided
    :return: Path to xyz file
    :rtype: Dict
    """

    assert scale_by in ['volume', 'a'], \
        'Only scaling by "a" and "volume" allowed'

    if  (scales is None and linspace is None and logspace is not None):
        scales = np.logspace(*logspace)
    elif (scales is None and linspace is not None and logspace is None):
        scales = np.linspace(*linspace)
    elif (scales is not None and linspace is None and logspace is None):
        pass
    else:
        raise ValueError(
            'Provide only one of factors, factor_logspacem factor_linspace'
        )

    atoms = get_atoms(**image)

    scales = np.array(scales)
    if scale_by == 'volume':
        scales = scales**(1/3)

    # Make a database of structures with the volumes scaled appropriately
    scaled_images = []
    for scale in scales:
        new_atoms = apply_scale(atoms, scale)
        scaled_images.append(new_atoms)

    scaled_images_file = 'scaled_unitcells_output.xyz'
    write(scaled_images_file, scaled_images, format='extxyz')

    return {'files': {'images': scaled_images_file}}
