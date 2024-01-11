#!/usr/bin/env python
'''
Relaxes the given atomic structure using ASE's built-in structure
optimizers
'''

from typing import Dict, Tuple, Optional
import numpy as np
import ase.optimize
from ase.io.trajectory import Trajectory
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms, get_logger

def atom_relax(
    calc_id: str,
    image: Dict,
    optimizer: str = 'GPMin', #GPMin is fast in many cases according to ASE docs
    properties: Tuple[str] = ('energy', 'forces'),
    fmax: float = 0.02,
    prefix: Optional[str] = None,
) -> Dict:
    """Relaxes the given tomic structure using ASE's built-in structure
    optimizers

    :param calc_id: calc_id specification
    :type calc_id: str
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param prefix: Prefix of output files, defaults to ''
    :type prefix: str, optional
    :param optimizer: ASE Optimizer class, defaults to 'GPMin'
    :type optimizer: str, optional
    :param fmax: Force convergence threshold in optimizer, defaults to 0.02
    :type fmax: float, optional
    :return: Dictionary of results
    :rtype: Dict
    """
    calc = load_calc(calc_id)
    atoms = get_atoms(**image)
    atoms.set_calculator(calc)
    logger = get_logger()

    if prefix is not None:
        prefix = prefix + '_'
    else:
        prefix = ''

    traj_file = prefix + 'atom_relax.traj'
    dyn = getattr(ase.optimize, optimizer)(atoms)
    traj = Trajectory(
        traj_file,
        'a',
        atoms,
        properties=properties
    )
    dyn.attach(traj)
    try:
        dyn.run(fmax=fmax)
    except Exception:
        logger.error('Failed to relax atoms')
        raise

    image_file = prefix + 'image_output.xyz'
    atoms.write(image_file, format='extxyz')

    energy = float(atoms.get_potential_energy())
    final_fmax = float(np.sqrt((atoms.get_forces() ** 2).sum(axis=1).max()))

    results = {
        'energy': energy,
        'final_fmax': final_fmax,
        'files':{
            'image': image_file,
            'traj': traj_file,
        }
    }
    return results
