#!/usr/bin/env python
'''
Calculates single point energy 

Author: mkphuthi@github.com
'''
import sys
from typing import Sequence
import numpy as np
import ase.optimize
from ase.io.trajectory import Trajectory
from asimtools.calculators import load_calc
from asimtools.job import Job
from asimtools.utils import (
    get_atoms,
    join_names,
    parse_command_line,
)

def atom_relax(
    calc_input: dict,
    image: dict = None,
    prefix: str = '',
    optimizer: str = 'GPMin', #GPMin is fastest according to docs
    properties: Sequence[str] = ('energy', 'forces'),
    fmax: float = 0.05,
    **kwargs
) -> dict:
    ''' 
    Relaxes the given structure using ASE
    '''
    sim_input = {'prefix': prefix}
    job = Job(calc_input, sim_input)
    job.start()

    calc = load_calc(calc_input)

    atoms = get_atoms(**image)
    atoms.set_calculator(calc)

    traj_file = join_names([prefix, 'output.traj'])
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
        job.update_status('failed')
        print('Failed to optimize atoms')
        raise

    image_file = join_names([prefix, 'image_output.xyz'])
    atoms.write(image_file, format='extxyz')

    energy = float(atoms.get_potential_energy())
    final_fmax = float(np.sqrt((atoms.get_forces() ** 2).sum(axis=1).max()))

    job.update_output({
        'energy': energy,
        'final_fmax': final_fmax
    })
    job.add_output_files({
        'image': image_file,
        'traj': traj_file,
    })

    job.complete()
    return job.get_output()


def main(argv):
    ''' main '''
    calc_input, sim_input = parse_command_line(argv)
    atom_relax(
        calc_input,
        **sim_input,
    )

    try:
        atom_relax(calc_input, **sim_input)
        sys.exit(0)
    except Exception:
        print('Failed to run singlepoint')
        raise


if __name__ == "__main__":
    main(sys.argv[1:])
