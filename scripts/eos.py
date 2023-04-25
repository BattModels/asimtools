#!/usr/bin/env python
'''
Calculates EOS

Author: mkphuthi@github.com
'''
from typing import TypeVar, Union, Tuple
import sys
from os.path import join
import os
import numpy as np
from asimtools.job import prepare_job, UnitJob, check_jobs
from asimtools.calculators import load_calc
from asimtools.utils import (
    write_yaml,
    write_csv_from_dict,
    read_yaml,
    get_atoms,
    join_names,
    parse_command_line,
)

Atoms = TypeVar('Atoms')


def eos(
    calc_params: dict,
    # atoms: Atoms = None,
    # symbol: str = None,
    # crystal_structure: str = None,
    # input_file: str = None,
    image: dict,
    prefix: str = '',
    nimages: int = 5,
    scale_range: Tuple[Union[float, float]] = (0.95, 1.05),
    # plot: bool = True,
    workdir: str = os.getcwd(),
    **kwargs
):
    '''
    Calculate EOS
    '''
    ########### Can do a helper function for this
    output_yaml = join_names([prefix, 'output.yaml'])
    output = {'status': 'started'}
    write_yaml(output_yaml, output)
    ###########

    # Fail early principle. There should usually be assert statements here
    # to make sure inputs are correct
    assert scale_range[0] < scale_range[1], 'x_scale[0] should be smaller than \
        x_scale[1]'

    atoms = get_atoms(**image)

    # Preprocessing the data, structures and data structures
    scales = np.linspace(scale_range[0], scale_range[1], nimages)

    # single point calculation script called here. Note that this will automatically 
    # submit nimages jobs
    sp_jobs = []
    for scale in scales:
        new_cell = atoms.get_cell() * scale
        new_atoms = atoms.copy()
        new_atoms.set_cell(new_cell)
        xprefix = join_names([prefix, f'x{scale:.2f}'])
        atoms_input_file = join_names([xprefix, 'input_atoms.xyz'])
        sp_sim_params = kwargs
        sp_sim_params.update({
            'script': '/Users/mphuthi/Documents/scrap/dev/asimtools/scripts/singlepoint.py',
            'image': {'input_file': atoms_input_file},
            'prefix': xprefix,
        })

        unitjob = UnitJob(
            calc_params=calc_params,
            sim_params=sp_sim_params,
            workdir=join(workdir, xprefix),
        )
        unitjob.gen_input_files()
        new_atoms.write(join(unitjob.get_workdir(), atoms_input_file))

        sp_jobs.append(unitjob)

    if kwargs.get('submit', True):
        for sp_job in sp_jobs:
            sp_job.submit()


    # Figure out if all the jobs are finished
    all_jobs_finished = check_jobs(sp_jobs)

    # If they are not finished, print the statues of each job which can be one 
    # of "In progress", "Complete" or "Failed"
    if not all_jobs_finished:
        sys.exit()

    # If all jobs are finished, prepare and write of results

    volumes = []
    energies = []
    for sp_job in sp_jobs:
        # You can get outputs directly from *output.yaml for generic outputs
        energy = sp_job.get_output().get('energy')
        energies.append(energy)

        # You can get outputs directly from *output_atoms.xyz for special cases
        volume = sp_job.get_output_images()[-1].get_volume()
        volumes.append(volume)


    results = {
        'scales': scales,
        'volumes': volumes,
        'energies': energies,
    }

    write_csv_from_dict(join_names([prefix, 'eos_output.csv']), results)

    ###### Can have a helper for this
    output.update({'status': 'complete'})
    write_yaml(output_yaml, output)
    results.update(output)
    ######
    return results


def main(argv):
    ''' Main '''
    calc_params, sim_params = parse_command_line(argv)
    eos(
        calc_params,
        **sim_params,
    )


if __name__ == "__main__":
    main(sys.argv[1:])
