#!/usr/bin/env python
'''
Relaxes structure using ASE while retaining symmetry. Requires ASE>=3.23

Author: mkphuthi@github.com
'''
from typing import Dict, Optional
from asimtools.asimmodules.geometry_optimization.symmetric_cell_relax import symmetric_cell_relax

def refined_symmetric_cell_relax(
    rough_calc_id: str,
    refined_calc_id: str,
    image: Dict,
    optimizer: str = 'BFGS',
    fmax: float = 0.003, #Roughly 0.48GPa
    optimizer_args: Optional[Dict] = None,
    fixsymmetry_args: Optional[Dict] = None,
    expcellfilter_args: Optional[Dict] = None,
) -> Dict:
    """Relaxes cell (and atoms) using ase.constraints.ExpCellFilter while retaining symmetry

    :param rough_calc_id: calc_id specification for cheap/rough calculator
    :type rough_calc_id: str
    :param refined_calc_id: calc_id specification for refining calculator
    :type refined_calc_id: str
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param optimizer: Any optimizer from ase.optimize, defaults to 'BFGS'
    :type optimizer: str, optional
    :param fmax: fmax in optimizer, defaults to 0.003 which is ~0.48GPa
    :type fmax: float, optional
    :param fixsymmetry_args: arguments for ase.constraints.FixSymmetry, defaults to None
    :type fixsymmetry_args: Optional[Dict], optional
    :param expcellfilter_args: arguments for ase.constraints.ExpCellFilter, defaults to None
    :type expcellfilter_args: Optional[Dict], optional
    :return: Results including relaxed structure energy
    :rtype: Dict
    """

    symmetric_cell_relax(
        calc_id=rough_calc_id,
        image=image,
        optimizer=optimizer,
        fmax=fmax,
        optimizer_args=optimizer_args,
        fixsymmetry_args=fixsymmetry_args,
        expcellfilter_args=expcellfilter_args,
    )

    results = symmetric_cell_relax(
        calc_id=refined_calc_id,
        image={'image_file': 'image_output.xyz'},
        optimizer=optimizer,
        fmax=fmax,
        optimizer_args=optimizer_args,
        fixsymmetry_args=fixsymmetry_args,
        expcellfilter_args=expcellfilter_args,
    )

    return results
