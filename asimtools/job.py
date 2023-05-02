'''
Set of tools for handling jobs in slurm, interactively or in the teminal for
a specific calculator.

Author: mkphuthi@github.com
'''

import subprocess
import os
from pathlib import Path
from ase.io import read
import numpy as np
from asimtools.utils import (
    read_yaml,
    write_yaml,
    join_names,
)


class Job():
    ''' Abstract class for the job object '''

    # pylint: disable=too-many-instance-attributes
    def __init__(
        self,
        calc_input: dict,
        sim_input: dict,
        workdir: str = None,
    ):

        self.calc_input = calc_input.copy()
        self.sim_input = sim_input.copy()
        if workdir is None:
            self.workdir = Path('.')
        else:
            self.workdir = Path(workdir)

    def mkworkdir(self):
        ''' Creates the work directory if it doesn't exist '''
        if not self.workdir.exists():
            self.workdir.mkdir()

    def start(self):
        ''' Updates the output to signal that the job was started '''
        self.update_output({'status': 'started'})

    def complete(self):
        ''' Updates the output to signal that the job was started '''
        self.update_output({'status': 'complete'})

    def update_status(self, status):
        ''' Updates job status to specificied value '''
        self.update_output({'status': status})

    def add_output_files(self, file_dict):
        ''' Adds a file to the output file list '''
        files = self.get_output().get('files', {})
        files.update(file_dict)
        self.update_output({'files': file_dict})

    def get_sim_input(self):
        ''' Get simulation input '''
        return self.sim_input

    def get_prefix(self):
        ''' Get current job prefix '''
        return self.sim_input.get('prefix', '')

    def get_output_yaml(self):
        ''' Get current output file '''
        prefix = self.get_prefix()
        out_fname = join_names([prefix, 'output.yaml'])
        output_yaml = self.workdir / out_fname
        return output_yaml

    def get_calc_input(self):
        ''' Get calculator input '''
        return self.calc_input

    def get_calc_params(self):
        ''' Get calculator parameters '''
        return self.calc_input.get('calc', {})

    def get_job_params(self):
        ''' Get job parameters '''
        return self.calc_input.get('job', {})

    def get_workdir(self):
        ''' Get working directory '''
        return self.workdir

    def get_output(self):
        ''' Get values in output.yaml '''
        output_yaml = self.get_output_yaml()
        if output_yaml.exists():
            return read_yaml(output_yaml)
        else:
            return {}

    def get_output_images(self, ext='xyz'):
        ''' Get all the output Atoms objects '''
        image_files = self.workdir / '*images_output.'+ext
        images = [read(image_file) for image_file in image_files]
        if len(images) > 0:
            return images

        return None

    def get_output_image(self, ext='extxyz'):
        ''' Get the output Atoms object '''
        files = self.get_output().get('files', {})
        image_file = files.get('image', False)
        if image_file:
            return read(self.workdir / image_file, format=ext)
        else:
            return None

    def set_sim_input(self, new_params):
        ''' Update simulation parameters '''
        self.sim_input.update(new_params)

    def set_calc_input(self, new_params):
        ''' Update calculator parameters '''
        self.calc_input.update(new_params)

    def set_workdir(self, workdir):
        ''' Set working directory '''
        self.workdir = workdir

    def check_job_status(self):
        ''' Check job status '''
        output = self.get_output()
        status = output.get('status', 'clean')
        complete = status == 'complete'
        return complete, status

    def update_output(self, output_update: dict):
        ''' Update output.yaml if it exists or write a new one '''
        output = self.get_output()
        output.update(output_update)
        write_yaml(self.get_output_yaml(), output)


class UnitJob(Job):
    ''' Unit job object with ability to submit slurm jobs '''
    def __init__(
        self,
        calc_input: dict,
        sim_input: dict,
        workdir = None,
    ):
        super().__init__(calc_input, sim_input, workdir)

        self.script = sim_input.get('script', None)

    def _gen_run_command(self):
        '''
        Generates the command to run the job, independent of slurm
        '''
        script = self.sim_input.get('script', False)
        assert script, 'Specify an existing script in sim_input'

        job_params = self.calc_input.get('job', {})
        run_prefix = job_params.get('run_prefix', '')
        run_suffix = job_params.get('run_suffix', '')

        txt = f'{run_prefix} '
        txt += f'{script} calc_input.yaml sim_input.yaml '
        txt += f'{run_suffix} '
        return txt

    def _gen_slurm_script(self, write: bool = True):
        '''
        Generates a slurm job file and if necessary writes it in the work
        directory for the job
        '''
        slurm_params = self.calc_input.get('slurm', {})

        txt = '#!/usr/bin/sh\n'
        for flag in slurm_params.get('flags'):
            txt += f'#SBATCH {flag}\n'
        txt += '\n'
        for command in slurm_params.get('precommands'):
            txt += command + '\n'
        txt += '\n'
        txt += 'echo "Job started on `hostname` at `date`"\n'
        txt += self._gen_run_command() + '\n'
        for command in slurm_params.get('postcommands'):
            txt += command + '\n'
        txt += 'echo "Job ended at `date`"'

        if write:
            slurm_file = self.workdir / 'job.sh'
            slurm_file.write_text(txt)

    def _gen_slurm_interactive_txt(self):
        '''
        Generates an srun command to run the job
        '''
        slurm_params = self.calc_input.get('slurm', False)

        txt = 'srun -u '
        txt += ' '.join(slurm_params.get('flags'))
        txt += self._gen_run_command()

        return txt

    def gen_input_files(self):
        ''' Write input files to working directory '''
        workdir = self.workdir
        sim_input = self.sim_input.copy()
        calc_input = self.calc_input
        job_params = calc_input.get('job', {})
        use_slurm = job_params.get('use_slurm', False)
        interactive = job_params.get('interactive', False)
        overwrite = self.sim_input.get('overwrite', False)
        if workdir.exists() and not overwrite:
            print(f'Skipped writing in {self.workdir}')
        else:
            if not workdir.exists():
                workdir.mkdir()
            write_yaml(workdir / 'sim_input.yaml', sim_input)
            write_yaml(workdir / 'calc_input.yaml', calc_input)
        if use_slurm and not interactive:
            self._gen_slurm_script()

    def submit(self):
        ''' Submit a job using slurm, interactively or in the terminal '''

        _, status = self.check_job_status()
        if status not in ('clean', 'discard'):
            return None

        if not self.sim_input.get('submit', False):
            return None

        cur_dir = Path('.').resolve()
        os.chdir(self.workdir)
        job_params = self.calc_input.get('job')
        sim_input = self.sim_input
        use_slurm = job_params.get('use_slurm', True)
        interactive = job_params.get('interactive', False)
        # Maybe want to do something clever with the exit code here
        if use_slurm and not interactive:
            command = ['sbatch', 'job.sh']
        elif use_slurm and interactive:
            command = self._gen_slurm_interactive_txt().split()
        else:
            command = self._gen_run_command().split()

        run_job = False
        overwrite = sim_input.get('overwrite', False)
        _, status = self.check_job_status()
        if overwrite or status == 'clean':
            run_job = True

        if run_job:
            completed_process = subprocess.run(
                command, check=False, capture_output=True, text=True
            )

            if completed_process.returncode != 0:
                err_txt = f'Failed to run {sim_input.get("script")} in '
                err_txt += f'{self.workdir} with command:\n'
                err_txt += f'{" ".join(command)}\n'
                err_txt += f'See {self.workdir / "errors.txt"} for details.'
                print(err_txt)
                with open('errors.txt', 'w', encoding='utf-8') as err_file:
                    err_file.write(completed_process.stderr)
                completed_process.check_returncode()

            #TODO: Get slurm job ID from stdout

        os.chdir(cur_dir)

        return None

def prepare_job(
    calc_input: dict,
    sim_input: dict,
    workdir: str,
):
    ''' Prepares a job's input files and gives a job object to handle '''
    job = UnitJob(calc_input, sim_input, workdir)
    job.gen_input_files()
    return job


def check_jobs(jobs: Job):
    ''' Checks status of jobs and prints a report '''

    completed = []
    print('+'*20)
    print('Job statuses: ')
    print('+'*20)
    for job in jobs:
        complete, status = job.check_job_status()
        print(f'{status}: {job.get_workdir()}')
        if complete or status == 'discard':
            completed.append(True)
        else:
            completed.append(False)

    if np.all(completed):
        return True
    else:
        return False
