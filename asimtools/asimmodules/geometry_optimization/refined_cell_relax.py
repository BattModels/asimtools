#!/usr/bin/env python
'''
Relaxes structure using ASE

Author: mkphuthi@github.com
'''

#pylint: disable=unused-argument
#pylint: disable=too-many-locals
#pylint: disable=too-many-arguments

from typing import Dict, Optional, Sequence
from asimtools.asimmodules.geometry_optimization.cell_relax import cell_relax

def refined_cell_relax(
    rough_calc_id: str,
    refined_calc_id: str,
    image: Dict,
    optimizer: str = 'BFGS',
    fmax: float = 0.002,
    mask: Optional[Sequence] = None,
    prefix: Optional[str] = None,
) -> Dict:
    """Relax cell using ASE Optimizer

    :param rough_calc_id: calc_id specification for cheap/rough calculator
    :type rough_calc_id: str
    :param refined_calc_id: calc_id specification for refining calculator
    :type refined_calc_id: str
    :param calc_id: calc_id specification
    :type calc_id: str
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param prefix: Prefix to output files, defaults to ''
    :type prefix: str, optional
    :param optimizer: Optimizer class to use from ase.optimize, defaults to 'BFGS'
    :type optimizer: str, optional
    :param fmax: Maximum stress in eV/$\AA^3$, defaults to 0.002
    :type fmax: float, optional
    :param mask: Mask to constrain cell deformation while relaxing, defaults to None
    :type mask: Optional[Sequence], optional
    :return: Dictionary of results including, final energy, stress and output files
    :rtype: Dict
    """

    cell_relax(
        calc_id=rough_calc_id,
        image=image,
        optimizer=optimizer,
        fmax=fmax,
        mask=mask,
        prefix=prefix,
    )

    results = cell_relax(
        calc_id=refined_calc_id,
        image={'image_file': 'image_output.xyz'},
        optimizer=optimizer,
        fmax=fmax,
        mask=mask,
        prefix=prefix,
    )

    return results
