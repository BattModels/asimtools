#!/usr/bin/env python
'''
Calculates single point energy 

Author: mkphuthi@github.com
'''
import sys
from typing import TypeVar
from asimtools.calculators import load_calc
from asimtools.job import Job
from asimtools.utils import (
    get_atoms,
    join_names,
    parse_command_line,
)

Atoms = TypeVar('Atoms')
# pylint: disable=too-many-arguments

def single_point(
    calc_input: dict,
    # atoms: Atoms = None,
    # symbol: str = None,
    # crystal_structure: str = None,
    # input_file: str = None,
    image: dict = None,
    prefix: str = '',
    **kwargs
):
    ''' 
    Calculates the single point energy, forces and stresses where possible
    '''
    sim_input = {'prefix': prefix}
    job = Job(calc_input, sim_input)
    job.start()

    calc = load_calc(calc_input)

    atoms = get_atoms(**image)
    atoms.set_calculator(calc)
    try:
        energy = atoms.get_potential_energy()
    except Exception:
        job.update_status('failed')
        raise

    image_file = join_names([prefix, 'image_output.xyz'])
    atoms.write(image_file, format='extxyz')

    job.update_output({'energy': float(energy)})
    job.add_output_files({'image': image_file})
    job.complete()
    return job.get_output()


def main(argv):
    ''' main '''
    calc_input, sim_input = parse_command_line(argv)
    single_point(
        calc_input,
        **sim_input,
    )

    try:
        single_point(calc_input, **sim_input)
        sys.exit(0)
    except Exception:
        print('Failed to run singlepoint')
        raise


if __name__ == "__main__":
    main(sys.argv[1:])
