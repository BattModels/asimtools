#!/usr/bin/env python
'''
Calculates EOS

Author: mkphuthi@github.com
'''

# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
from typing import Tuple, Dict, List
import numpy as np
from ase import Atoms
from asimtools.job import UnitJob
from asimtools.utils import (
    get_atoms,
    join_names,
    read_yaml,
)


def singlepoint_jobs(
    atoms: Atoms,
    scales: List,
    prefix: str,
    calc_input: Dict,
    **kwargs
) -> List[UnitJob]:
    """
    Generate list of singlepoint calculations
    """
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

    return sp_jobs


def eos(
    calc_input: dict,
    image: dict,
    prefix: str = '',
    nimages: int = 5,
    scale_range: Tuple[float] = (0.95, 1.05),
    # plot: bool = True,
    # workdir: str = '.',
    **kwargs
) -> Dict:
    '''
    Calculate EOS
    '''

    # Fail early principle. There should usually be assert statements here
    # to make sure inputs are correct
    assert scale_range[0] < scale_range[1], 'x_scale[0] should be smaller than\
         x_scale[1]'

    atoms = get_atoms(**image)

    # Preprocessing the data, structures and data structures
    scales = np.linspace(scale_range[0], scale_range[1], nimages)

    # single point calculation script called here. Note that this will
    # automatically submit nimages jobs

    calculator_inputs = calc_input.get('calculator_filenames')
    if isinstance(calculator_inputs,list):
        eos_input = read_yaml(calculator_inputs.pop(0))
    else:
        eos_input = calc_input

    sp_jobs = singlepoint_jobs(
        atoms,
        scales,
        prefix,
        eos_input,
        **kwargs
    )
    work_directories = [str(sp_job.workdir.resolve()) for sp_job in sp_jobs]
    
    eos_postprocess_job_input = {
        **kwargs,
        'script': 'eos_postprocess.py',
        'prefix': 'postprocess',
        'sp_work_directories': work_directories,
        'scales': scales.tolist(),
        # **kwargs
    }

    if isinstance(calculator_inputs,list):
        postprocess_input = read_yaml(calculator_inputs.pop(0))
    else:
        postprocess_input = calc_input
    
    eos_postprocess_job = UnitJob(
        postprocess_input,
        eos_postprocess_job_input,
        eos_postprocess_job_input['prefix']
    )
    eos_postprocess_job.gen_input_files()

    if kwargs.get('submit', True):
        job_ids = []
        for sp_job in sp_jobs:
            job_ids.append(sp_job.submit())
        print(job_ids)
        eos_postprocess_job.submit(dependency=job_ids)
