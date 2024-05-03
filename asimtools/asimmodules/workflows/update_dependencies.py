#!/usr/bin/env python
'''
Mostly internal tool used to chain together batch jobs that have nested
dependencies

Author: mkphuthi@github.com
'''

from typing import Dict, Optional
from pathlib import Path
import os
import subprocess
import logging
from asimtools.utils import read_yaml

def update_dependencies(
    prev_step_dir: Optional[os.PathLike] = None,
    next_step_dir: Optional[os.PathLike] = None,
    skip_failed: bool = False,
) -> Dict:
    """Asimmodule for connecting jobs where an asimmodule launches slurm jobs
    internally, examples are image_array, sim_array and calc_array

    :param prev_step_dir: workdir of previous step, usually step-[N-1] where N
        is current step, defaults to None
    :type prev_step_dir: Optional[os.PathLike], optional
    :param next_step_dir: workdir of next step, usually step-[N+1] where N
        is current step, defaults to None
    :type next_step_dir: Optional[os.PathLike], optional
    :param skip_failed: whether to submit if any jobs failed or not,
        defaults to False
    :type skip_failed: bool, optional
    :return: Nothing
    :rtype: Dict
    """

    prev_step_dir = Path(prev_step_dir)
    next_step_dir = Path(next_step_dir)
    prev_output = read_yaml(prev_step_dir / 'output.yaml')
    next_output = read_yaml(next_step_dir / 'output.yaml')
    next_job_ids = next_output.get('job_ids', [])
    job_ids = prev_output.get('job_ids', None)
    if job_ids is not None:
        job_ids = [str(job_id) for job_id in job_ids]
        start_cond = 'afterok'
        if skip_failed:
            start_cond = 'afterany'

        for next_job_id in next_job_ids:
            command = f'scontrol update JobId={next_job_id} '
            command += f'dependency={start_cond}:{":".join(job_ids)}'
            command = command.split()

            completed_process = subprocess.run(
                command, check=False, capture_output=True, text=True,
            )

            with open('sbatch_stdout.txt', 'w', encoding='utf-8') as f:
                f.write(completed_process.stdout)

            if completed_process.returncode != 0:
                err_txt = f'Failed to run {" ".join(command)}\n'
                err_txt += 'See sbatch.stderr.txt for details.'
                logging.error(err_txt)
                with open('sbatch_stderr.txt', 'w', encoding='utf-8') as f:
                    f.write(completed_process.stderr)
                completed_process.check_returncode()
            return {}

    return {}
