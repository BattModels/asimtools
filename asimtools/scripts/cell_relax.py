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
    prefix: str = '',
    optimizer: str = 'BFGS',
    smax: float = 0.002,
    mask: Optional[Sequence] = None,
) -> Dict:
    ''' 
    Relaxes the given lattice using ASE up to the threshold smax,
    needs stresses to be implemented in calculator.
    '''
    calc = load_calc(calc_id)
    atoms = get_atoms(**image)
    atoms.set_calculator(calc)

    traj_file = join_names([prefix, 'cell_relax.traj'])
    sf = StrainFilter(atoms, mask=mask)
    dyn = getattr(ase.optimize, optimizer)(sf)
    traj = Trajectory(
        traj_file,
        'w',
        atoms,
        properties=['energy', 'forces', 'stress'],
    )
    dyn.attach(traj)
    try:
        dyn.run(fmax=smax)
    except Exception:
        print('Failed to optimize atoms')
        raise

    image_file = join_names([prefix, 'image_output.xyz'])
    atoms.write(image_file, format='extxyz')

    energy = float(atoms.get_potential_energy())
    final_smax = float(atoms.get_stress().max())

    results = {
        'energy': energy,
        'final_smax': final_smax,
        'files':{
            'image': image_file,
            'traj': traj_file,
        }
    }
    return results
