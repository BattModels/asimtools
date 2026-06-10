#!/usr/bin/env python
'''
Relaxes both atoms and cell coordinates simultaneously.

Author: mkphuthi@github.com
'''
from typing import Dict, Optional
import ase.optimize
from ase.filters import ExpCellFilter
from ase.io.trajectory import Trajectory
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms, write_atoms, get_logger, join_names

def optimize(
    calculator: Dict,
    image: Dict,
    optimizer: str = 'LBFGS',
    fmax: float = 0.003, #Roughly 0.48GPa
    optimizer_args: Optional[Dict] = None,
    steps: int = 1000,
    expcellfilter_args: Optional[Dict] = None,
    prefix: Optional[str] = None,
) -> Dict:
    """Relaxes cell (and atoms) using ase.constraints.ExpCellFilter while retaining symmetry

    :param calculator: Calculator specification, see :func:`asimtools.calculators.load_calc`
    :type calculator: Dict
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param optimizer: Any optimizer from ase.optimize, defaults to 'BFGS'
    :type optimizer: str, optional
    :param fmax: fmax in optimizer, defaults to 0.003 which is ~0.48GPa
    :type fmax: float, optional
    :param optimizer_args: arguments for optimizer, defaults to None
    :type optimizer_args: Optional[Dict], optional
    :param steps: Max. number of steps to run, defaults to 1000
    :type steps: int, optional
    :param expcellfilter_args: arguments for ase.constraints.ExpCellFilter, defaults to None
    :type expcellfilter_args: Optional[Dict], optional
    :return: Results including relaxed structure energy
    :rtype: Dict
    """
    logger = get_logger()
    if expcellfilter_args is None:
        expcellfilter_args = {}
    if optimizer_args is None:
        optimizer_args = {}

    if prefix is not None:
        prefix = prefix + '_'
    else:
        prefix = ''

    calc = load_calc(calculator=calculator)
    atoms = get_atoms(**image)
    atoms.calc = calc

    ecf = ExpCellFilter(atoms, **expcellfilter_args)

    dyn = getattr(ase.optimize, optimizer)(ecf)
    traj_file = join_names([prefix, 'optimize.traj'])[:-2] # Don't include __
    traj = Trajectory(
        traj_file,
        'a', # If there are multiple relaxations, we keep all trajectories
        atoms,
        properties=['energy', 'forces', 'stress'],
    )
    dyn.attach(traj)

    try:
        converged = dyn.run(fmax=fmax, steps=steps)
    except Exception:
        print('Failed to optimize atoms')
        raise

    energy = float(atoms.get_potential_energy())
    final_fmax = float(np.sqrt((atoms.get_forces() ** 2).sum(axis=1).max()))
    final_smax = float(atoms.get_stress().max())

    if converged:
        image_file = prefix + 'image_output.xyz'
        write_atoms(
            image_file,
            atoms,
            format='extxyz',
        )

        results = {
            'energy': energy,
            'final_fmax': final_fmax,
            'final_smax': final_smax,
            'converged': converged,
            'files':{
                'image': image_file,
                'traj': traj_file,
            }
        }
        return results
    else:
        logger.error('Failed to relax atoms below fmax={}'.format(final_fmax))
        raise ValueError(f'Only relaxed below fmax={final_fmax}')