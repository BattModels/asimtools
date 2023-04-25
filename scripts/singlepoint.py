#!/usr/bin/env python
'''
Calculates single point energy 

Author: mkphuthi@github.com
'''
import sys
from typing import TypeVar
from asimtools.job import prepare_job
from asimtools.calculators import load_calc
from asimtools.utils import (
    write_yaml,
    read_yaml,
    get_atoms,
    join_names,
    parse_command_line,
)

Atoms = TypeVar('Atoms')
# pylint: disable=too-many-arguments

def single_point(
    calc_params: dict,
    # atoms: Atoms = None,
    # symbol: str = None,
    # crystal_structure: str = None,
    # input_file: str = None,
    image: dict = None,
    prefix: str = '',
    **kwargs,
):
    ''' 
    Calculates the single point energy, forces and stresses where possible
    '''
    output_yaml = join_names([prefix, 'output.yaml'])

    calc = load_calc(calc_params)
    atoms = get_atoms(**image)
    atoms.set_calculator(calc)
    results = {'status': 'started'}
    write_yaml(output_yaml, results)
    try:
        energy = atoms.get_potential_energy()
    except Exception:
        results['status'] = 'failed'
        write_yaml(output_yaml, results)
        raise

    atoms.write(join_names([prefix, 'output_atoms.xyz']), format='extxyz')
    results.update({
        'energy': float(energy),
        'status': 'complete',
    })
    write_yaml(output_yaml, results)

    return results


def main(argv):
    ''' main '''
    calc_params, sim_params = parse_command_line(argv)
    single_point(
        calc_params,
        **sim_params,
    )

    try:
        single_point(calc_params, **sim_params)
        sys.exit(0)
    except Exception:
        print('Failed to run singlepoint')
        raise


if __name__ == "__main__":
    main(sys.argv[1:])
