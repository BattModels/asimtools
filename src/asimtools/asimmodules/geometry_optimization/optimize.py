#!/usr/bin/env python
'''
Relaxes both atomic positions and cell simultaneously using ASE's
FrechetCellFilter. Use symmetric_cell_relax to additionally retain symmetry.

Author: mkphuthi@github.com
'''
from typing import Dict, Optional
from ase.filters import FrechetCellFilter
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms, get_logger
from asimtools.asimmodules.geometry_optimization.utils import run_optimization

def optimize(
    calculator: Dict,
    image: Dict,
    optimizer: str = 'LBFGS',
    fmax: float = 0.003, #Roughly 0.48GPa
    optimizer_args: Optional[Dict] = None,
    steps: int = 1000,
    cellfilter_args: Optional[Dict] = None,
    prefix: Optional[str] = None,
) -> Dict:
    """Relaxes cell and atomic positions using ase.filters.FrechetCellFilter

    :param calculator: Calculator specification, see :func:`asimtools.calculators.load_calc`
    :type calculator: Dict
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param optimizer: Any optimizer class name from ase.optimize, defaults to 'LBFGS'
    :type optimizer: str, optional
    :param fmax: fmax in optimizer, defaults to 0.003 which is ~0.48GPa
    :type fmax: float, optional
    :param optimizer_args: arguments for optimizer, defaults to None
    :type optimizer_args: Optional[Dict], optional
    :param steps: Max. number of steps to run, defaults to 1000
    :type steps: int, optional
    :param cellfilter_args: arguments for ase.filters.FrechetCellFilter, defaults to None
    :type cellfilter_args: Optional[Dict], optional
    :param prefix: Prefix of output files, defaults to None
    :type prefix: str, optional
    :return: Results including relaxed structure energy, fmax, stress and files
    :rtype: Dict
    """
    if cellfilter_args is None:
        cellfilter_args = {}

    logger = get_logger()
    calc = load_calc(calculator=calculator)
    atoms = get_atoms(**image)
    atoms.calc = calc

    ecf = FrechetCellFilter(atoms, **cellfilter_args)

    return run_optimization(
        atoms,
        optimizer=optimizer,
        fmax=fmax,
        steps=steps,
        optimize_obj=ecf,
        prefix=prefix,
        label='optimize',
        properties=('energy', 'forces', 'stress'),
        optimizer_args=optimizer_args,
        report_stress=True,
        logger=logger,
    )
