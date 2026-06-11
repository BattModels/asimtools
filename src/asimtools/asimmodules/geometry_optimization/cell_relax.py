#!/usr/bin/env python
'''
Relaxes the cell of a structure using ASE's StrainFilter, keeping the
fractional atomic coordinates fixed.

Author: mkphuthi@github.com
'''

from typing import Dict, Optional, Sequence
from ase.filters import StrainFilter
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms, get_logger
from asimtools.asimmodules.geometry_optimization.utils import run_optimization

def cell_relax(
    calculator: Dict,
    image: Dict,
    optimizer: str = 'LBFGS',
    fmax: float = 0.002,
    steps: int = 1000,
    mask: Optional[Sequence] = None,
    optimizer_args: Optional[Dict] = None,
    prefix: Optional[str] = None,
) -> Dict:
    """Relax cell using ASE's StrainFilter

    :param calculator: Calculator specification, see :func:`asimtools.calculators.load_calc`
    :type calculator: Dict
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param optimizer: Optimizer class name from ase.optimize, defaults to 'LBFGS'
    :type optimizer: str, optional
    :param fmax: Maximum stress in eV/$\\AA^3$, defaults to 0.002
    :type fmax: float, optional
    :param steps: Maximum number of steps to run, defaults to 1000
    :type steps: int, optional
    :param mask: Mask to constrain cell deformation while relaxing, defaults to None
    :type mask: Optional[Sequence], optional
    :param optimizer_args: Extra arguments for the optimizer class, defaults to None
    :type optimizer_args: Optional[Dict], optional
    :param prefix: Prefix to output files, defaults to None
    :type prefix: str, optional
    :return: Dictionary of results including final energy, stress and output files
    :rtype: Dict
    """
    logger = get_logger()
    calc = load_calc(calculator=calculator)
    atoms = get_atoms(**image)
    atoms.calc = calc

    sf = StrainFilter(atoms, mask=mask)

    return run_optimization(
        atoms,
        optimizer=optimizer,
        fmax=fmax,
        steps=steps,
        optimize_obj=sf,
        prefix=prefix,
        label='cell_relax',
        properties=('energy', 'forces', 'stress'),
        optimizer_args=optimizer_args,
        report_stress=True,
        logger=logger,
    )
