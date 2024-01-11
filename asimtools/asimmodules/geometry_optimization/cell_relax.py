#!/usr/bin/env python
'''
Relaxes structure using ASE

Author: mkphuthi@github.com
'''

#pylint: disable=unused-argument
#pylint: disable=too-many-locals
#pylint: disable=too-many-arguments

from typing import Dict, Optional, Sequence
import ase.optimize
from ase.constraints import StrainFilter
from ase.io.trajectory import Trajectory
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms, join_names

def cell_relax(
    calc_id: str,
    image: Dict,
    optimizer: str = 'BFGS',
    fmax: float = 0.002,
    mask: Optional[Sequence] = None,
    prefix: Optional[str] = None,
) -> Dict:
    """Relax cell using ASE Optimizer

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
    calc = load_calc(calc_id)
    atoms = get_atoms(**image)
    atoms.set_calculator(calc)

    if prefix is not None:
        prefix = prefix + '_'
    else:
        prefix = ''

    traj_file = join_names([prefix, 'cell_relax.traj'])[:-2] # Don't include __
    sf = StrainFilter(atoms, mask=mask)
    dyn = getattr(ase.optimize, optimizer)(sf)
    traj = Trajectory(
        traj_file,
        'a',
        atoms,
        properties=['energy', 'forces', 'stress'],
    )
    dyn.attach(traj)
    try:
        dyn.run(fmax=fmax)
    except Exception:
        print('Failed to optimize atoms')
        raise

    image_file = join_names([prefix, 'image_output.xyz'])[:-2]
    atoms.write(image_file, format='extxyz')

    energy = float(atoms.get_potential_energy())
    final_fmax = float(atoms.get_stress().max())

    results = {
        'energy': energy,
        'final_fmax': final_fmax,
        'files':{
            'image': image_file,
            'traj': traj_file,
        }
    }
    return results
