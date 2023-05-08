#!/usr/bin/env python
'''
Calculates EOS

Author: mkphuthi@github.com
'''

#pylint: disable=too-many-arguments
from typing import TypeVar
import sys
from asimtools.job import UnitJob, Job, check_jobs
from asimtools.utils import (
    get_images,
    parse_command_line,
)

Atoms = TypeVar('Atoms')


def multiple_images(
    calc_input: dict,
    images: dict,
    subscript_input: dict,
    prefix: str = '',
    workdir: str = '.',
    **kwargs
):
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


def main(argv):
    ''' Main '''
    calc_input, sim_input = parse_command_line(argv)
    multiple_images(
        calc_input,
        **sim_input,
    )


if __name__ == "__main__":
    main(sys.argv[1:])
