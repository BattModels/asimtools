#!/usr/bin/env python
'''
Calculates EOS

Author: mkphuthi@github.com
'''

#pylint: disable=too-many-arguments
#pylint: disable=too-many-locals
from typing import TypeVar, Tuple, Dict
import sys
import numpy as np
from asimtools.job import UnitJob, Job, check_jobs
from asimtools.utils import (
    write_csv_from_dict,
    get_atoms,
    join_names,
)

Atoms = TypeVar('Atoms')


def eos(
    calc_input: dict,
    image: dict,
    prefix: str = '',
    nimages: int = 5,
    scale_range: Tuple[float] = (0.95, 1.05),
    # plot: bool = True,
    workdir: str = '.',
    **kwargs
) -> Dict:
    '''
    Calculate EOS
    '''
    job = Job(calc_input, {'prefix': prefix, 'workdir': workdir})
    job.start()

    # Fail early principle. There should usually be assert statements here
    # to make sure inputs are correct
    assert scale_range[0] < scale_range[1], 'x_scale[0] should be smaller than\
         x_scale[1]'

    atoms = get_atoms(**image)

    # Preprocessing the data, structures and data structures
    scales = np.linspace(scale_range[0], scale_range[1], nimages)

    # single point calculation script called here. Note that this will
    # automatically submit nimages jobs
    sp_jobs = []
    for scale in scales:
        new_cell = atoms.get_cell() * scale
        new_atoms = atoms.copy()
        new_atoms.set_cell(new_cell)
        xprefix = join_names([prefix, f'x{scale:.2f}'])
        sp_sim_input = kwargs
        sp_sim_input.update({
            'script': 'singlepoint.py',
            # 'script': 'asim-run',
            'prefix': xprefix,
        })

        unitjob = UnitJob(
            calc_input,
            sp_sim_input,
            xprefix,
        )
        unitjob.set_input_image(new_atoms)
        unitjob.gen_input_files()
        sp_jobs.append(unitjob)

    if kwargs.get('submit', True):
        for sp_job in sp_jobs:
            sp_job.submit()

    # Figure out if all the jobs are finished
    all_jobs_finished = check_jobs(sp_jobs)

    # If they are not finished, print the statues of each job
    if not all_jobs_finished:
        sys.exit()

    # If all jobs are finished, prepare and write of results
    volumes = []
    energies = []
    for sp_job in sp_jobs:
        # You can get outputs directly from *output.yaml for saved outputs
        energy = sp_job.get_output().get('energy')
        energies.append(energy)

        # You can get outputs directly from *image_output.xyz for special cases
        volume = sp_job.get_output_image().get_volume()
        volumes.append(volume)

    results = {
        'scales': scales,
        'volumes': volumes,
        'energies': energies,
    }

    eos_table_file = join_names([prefix, 'eos_output.csv'])
    write_csv_from_dict(eos_table_file, results)

    job.add_output_files({'eos': eos_table_file})
    job.complete()
    results.update(job.get_output())
    return results
