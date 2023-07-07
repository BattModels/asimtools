#!/usr/bin/env python
'''
Relaxes structure using ASE

Author: mkphuthi@github.com
'''

#pylint: disable=unused-argument
#pylint: disable=too-many-locals
#pylint: disable=too-many-arguments

from typing import Dict, Tuple
import numpy as np
import ase.optimize
from ase.io.trajectory import Trajectory
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms, join_names

def atom_relax(
    calc_id: str,
    image: Dict,
    optimizer: str = 'GPMin', #GPMin is fastest according to ASE docs
    properties: Tuple[str] = ('energy', 'forces'),
    fmax: float = 0.02,
) -> Dict:
    ''' 
    Relaxes the given structure using ASE
    '''
    calc = load_calc(calc_id)
    atoms = get_atoms(**image)
    atoms.set_calculator(calc)

    traj_file = 'atom_relax.traj'
    dyn = getattr(ase.optimize, optimizer)(atoms)
    traj = Trajectory(
        traj_file,
        'w',
        atoms,
        properties=properties
    )
    dyn.attach(traj)
    try:
        dyn.run(fmax=fmax)
    except Exception:
        print('Failed to optimize atoms')
        raise

    image_file = 'image_output.xyz'
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
