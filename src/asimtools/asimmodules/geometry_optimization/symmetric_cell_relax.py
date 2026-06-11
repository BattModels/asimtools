#!/usr/bin/env python
'''
Relaxes cell and atomic positions using ASE while retaining symmetry.
Requires ASE>=3.23

Author: mkphuthi@github.com
'''
from typing import Dict, Optional
from ase.filters import FrechetCellFilter
from ase.constraints import FixSymmetry
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms, get_logger
from asimtools.asimmodules.geometry_optimization.utils import run_optimization

def symmetric_cell_relax(
    calculator: Dict,
    image: Dict,
    optimizer: str = 'LBFGS',
    fmax: float = 0.003, #Roughly 0.48GPa
    steps: int = 1000,
    optimizer_args: Optional[Dict] = None,
    fixsymmetry_args: Optional[Dict] = None,
    cellfilter_args: Optional[Dict] = None,
    prefix: Optional[str] = None,
) -> Dict:
    """Relaxes cell and atomic positions using ase.filters.FrechetCellFilter
    while retaining symmetry

    :param calculator: Calculator specification, see :func:`asimtools.calculators.load_calc`
    :type calculator: Dict
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param optimizer: Any optimizer class name from ase.optimize, defaults to 'LBFGS'
    :type optimizer: str, optional
    :param fmax: fmax in optimizer, defaults to 0.003 which is ~0.48GPa
    :type fmax: float, optional
    :param steps: Max. number of steps to run, defaults to 1000
    :type steps: int, optional
    :param optimizer_args: arguments for optimizer, defaults to None
    :type optimizer_args: Optional[Dict], optional
    :param fixsymmetry_args: arguments for ase.constraints.FixSymmetry, defaults to None
    :type fixsymmetry_args: Optional[Dict], optional
    :param cellfilter_args: arguments for ase.filters.FrechetCellFilter, defaults to None
    :type cellfilter_args: Optional[Dict], optional
    :param prefix: Prefix of output files, defaults to None
    :type prefix: str, optional
    :return: Results including relaxed structure energy, fmax, stress and files
    :rtype: Dict
    """
    if fixsymmetry_args is None:
        fixsymmetry_args = {}
    if cellfilter_args is None:
        cellfilter_args = {}

    logger = get_logger()
    calc = load_calc(calculator=calculator)
    atoms = get_atoms(**image)
    atoms.calc = calc

    atoms.set_constraint(FixSymmetry(atoms, **fixsymmetry_args))
    ecf = FrechetCellFilter(atoms, **cellfilter_args)

    return run_optimization(
        atoms,
        optimizer=optimizer,
        fmax=fmax,
        steps=steps,
        optimize_obj=ecf,
        prefix=prefix,
        label='symmetric_cell_relax',
        properties=('energy', 'forces', 'stress'),
        optimizer_args=optimizer_args,
        report_stress=True,
        logger=logger,
    )
