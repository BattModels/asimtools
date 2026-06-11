#!/usr/bin/env python
'''
Relaxes the given atomic structure using ASE's built-in structure
optimizers

Author: mkphuthi@github.com
'''

from typing import Dict, Sequence, Optional
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms, get_logger
from asimtools.asimmodules.geometry_optimization.utils import run_optimization

def atom_relax(
    calculator: Dict,
    image: Dict,
    optimizer: str = 'LBFGS',
    properties: Sequence[str] = ('energy', 'forces'),
    fmax: float = 0.02,
    steps: int = 1000,
    optimizer_args: Optional[Dict] = None,
    prefix: Optional[str] = None,
) -> Dict:
    """Relaxes the given atomic structure using ASE's built-in structure
    optimizers

    :param calculator: Calculator specification, see :func:`asimtools.calculators.load_calc`
    :type calculator: Dict
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param optimizer: ASE Optimizer class name from ase.optimize, defaults to 'LBFGS'
    :type optimizer: str, optional
    :param properties: Properties stored in the trajectory, defaults to ('energy', 'forces')
    :type properties: Sequence[str], optional
    :param fmax: Force convergence threshold in optimizer, defaults to 0.02
    :type fmax: float, optional
    :param steps: Maximum number of steps to run, defaults to 1000
    :type steps: int, optional
    :param optimizer_args: Extra arguments for the optimizer class, defaults to None
    :type optimizer_args: Optional[Dict], optional
    :param prefix: Prefix of output files, defaults to None
    :type prefix: str, optional
    :return: Dictionary of results
    :rtype: Dict
    """
    logger = get_logger()
    calc = load_calc(calculator=calculator)
    atoms = get_atoms(**image)
    atoms.calc = calc

    return run_optimization(
        atoms,
        optimizer=optimizer,
        fmax=fmax,
        steps=steps,
        prefix=prefix,
        label='atom_relax',
        properties=properties,
        optimizer_args=optimizer_args,
        logger=logger,
    )
