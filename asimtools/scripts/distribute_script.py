#!/usr/bin/env python
'''
Apply the same script on multiple images/calc inputs

Author: mkphuthi@github.com
'''

from typing import TypeVar, Dict
import sys
from asimtools.job import UnitJob, Job, check_jobs
from asimtools.utils import get_images


def distribute_script(
    calc_input: Dict,
    script: str,
    images: Dict = None,
    image: Dict = None,
    prefix: str = '',
    workdir: str = '.',
    **kwargs
) -> Dict:
    '''
    Repeat a script on multiple images
    '''
    job = Job(calc_input, {'prefix': prefix, 'workdir': workdir})
    job.start()

    # Fail early principle. There should usually be assert statements here
    # to make sure inputs are correct

    if images is not None:
        images = get_images(**images)
    else:
        images = [get_atoms(**image)]

    configs = 
    if 
    unitjobs = []
    for ind, atoms in enumerate(images):
        unitjob = UnitJob(
            calc_input,
            subscript_input,
            workdir=f'id-{ind}',
        )
        unitjob.set_input_image(atoms)
        unitjobs.append(unitjob)

    for unitjob in unitjobs:
        if kwargs.get('submit', True):
            unitjob.submit()
        else:
            unitjob.gen_input_files()
            print(unitjob.get_workdir().resolve())

    # Figure out if all the jobs are finished
    all_jobs_finished = check_jobs(unitjobs)

    # If they are not finished, print the statues of each job
    if not all_jobs_finished:
        # Can do clever things with error codes here potentially
        sys.exit(1)

    job.complete()

    return job.get_output()
