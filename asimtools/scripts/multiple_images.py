#!/usr/bin/env python
'''
Apply the same script on multiple images

Author: mkphuthi@github.com

This is a breaking case of the current architecture, will need to 
rethink fundamental parts of the Job object and how it handles jobs
with different submission requirements e.g. a job needing one core
vs 32 cores. Basically a way to hand scripts calling other scripts
at various levels.
'''

from typing import TypeVar, Dict
import sys
from asimtools.job import UnitJob, Job, check_jobs
from asimtools.utils import get_images

Atoms = TypeVar('Atoms')


def multiple_images(
    calc_input: Dict,
    images: Dict,
    subscript_input: Dict,
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

    images = get_images(**images)

    unitjobs = []
    for ind, atoms in enumerate(images):
        unitjob = UnitJob(
            calc_input,
            subscript_input,
            workdir=f'id_{ind}',
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
