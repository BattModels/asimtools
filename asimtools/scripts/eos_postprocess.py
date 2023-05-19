#!/usr/bin/env python
'''
Calculates EOS

Author: mkphuthi@github.com
'''

#pylint: disable=too-many-arguments
#pylint: disable=too-many-locals
from typing import Dict, List
import sys

from asimtools.job import  Job, check_jobs, load_jobs_from_directory
from asimtools.utils import (
    write_csv_from_dict,
    get_atoms,
    join_names,
)

def eos_postprocess(
    calc_input: dict,
    scales: list,
    sp_work_directories: list,
    prefix: str = '',
    # plot: bool = True,
    workdir: str = '',
    **kwargs
) -> Dict:
    '''
    Calculate EOS after singlepoint calculations are completed
    '''
    job = Job(calc_input, {'prefix': prefix, 'workdir': workdir})
    job.start()

    
    sp_jobs = load_jobs_from_directory(
        sp_work_directories
    )
    print(sp_jobs)
    all_jobs_finished = check_jobs(sp_jobs)
    print(all_jobs_finished)
    # If they are not finished, print the statues of each job
    # if not all_jobs_finished:
    #     sys.exit()

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
    print(results)
    eos_table_file = join_names([prefix, 'eos_output.csv'])
    print(eos_table_file)
    write_csv_from_dict(eos_table_file, results)

    job.add_output_files({'eos': eos_table_file})
    job.complete()
    results.update(job.get_output())
    return results
