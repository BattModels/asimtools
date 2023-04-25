'''
Set of tools for handling jobs in slurm, interactively or in the teminal for 
a specific calculator.

Author: mkphuthi@github.com
'''

import subprocess
from abc import ABC
import os
from os.path import join
from glob import glob
from ase.io import read
import numpy as np
from asimtools.utils import (
    read_yaml,
    write_yaml,
    join_names,
)

class Job(ABC):
    ''' Abstract class for the job object '''

    # pylint: disable=too-many-instance-attributes
    def __init__(
        self,
        sim_params: dict,
        calc_params: dict,
        workdir: str,
    ):

        self.calc_params = calc_params
        self.slurm = calc_params.get('slurm', True)
        self.interactive = calc_params.get('interactive', False)
        self.batch = calc_params.get('use_srun', False)
        if self.slurm:
            assert (self.interactive and not self.batch) or \
                (not self.interactive and self.batch), \
                     'If using slurm, specify one of batch or interactive job'

        self.sim_params = sim_params
        self.workdir = workdir
        self.script = sim_params.get('script', None)
        self.overwrite = sim_params.get('overwrite', False)
        prefix = self.sim_params.get('prefix', '')
        out_fname = join_names([prefix, 'output.yaml'])
        self.output_file = join(self.workdir, out_fname)

        self.ready = False   # Ready to submit?

    def mkworkdir(self):
        ''' Creates the work directory if it doesn't exist '''
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

    def update_output(self, d: dict):
        ''' Update output.yaml '''
        output = read_yaml(self.output_file)
        output.update(d)
        write_yaml(self.output_file, output)

    def get_sim_params(self):
        ''' Get simulation parameters '''
        return self.sim_params

    def get_calc_params(self):
        ''' Get calculator parameters '''
        return self.calc_params

    def get_workdir(self):
        ''' Get working directory '''
        return self.workdir

    def get_output(self):
        ''' Get values in output.yaml '''
        return read_yaml(self.output_file)

    def get_output_images(self, ext='xyz'):
        ''' Get the output Atoms object '''
        image_files = glob(join(self.workdir, '*output_atoms.'+ext))
        images = [read(image_file) for image_file in image_files]
        return images

    def set_sim_params(self, new_params):
        ''' Update simulation parameters '''
        self.sim_params.update(new_params)

    def set_calc_params(self, new_params):
        ''' Update calculator parameters '''
        self.calc_params.update(new_params)

    def set_workdir(self, workdir):
        ''' Set working directory '''
        self.workdir = workdir


class UnitJob(Job):
    ''' Unit job object with ability to submit slurm jobs '''
    def __init__(
        self,
        sim_params: dict,
        calc_params: dict,
        workdir: str,
    ):
        super(UnitJob, self).__init__(sim_params, calc_params, workdir)

        self.run_cmd = self._gen_run_command()
        if self.calc_params.get('slurm', False):
            if self.calc_params.get('interactive', False):
                self.interactive_cmd = self._gen_slurm_interactive_txt()

    def _gen_run_command(self):
        '''
        Generates the command to run the job, independent of slurm
        '''
        script = self.sim_params.get('script', False)
        assert script, 'Specify an existing script in sim_params'

        run_prefix = self.calc_params.get('run_prefix', '')
        run_suffix = self.calc_params.get('run_suffix', '')

        txt = f'{run_prefix} '
        txt += f'python {script} calc_params.yaml sim_params.yaml '
        txt += f'{run_suffix} '
        return txt

    def _gen_slurm_script(self, write: bool = True):
        '''
        Generates a slurm job file and if necessary writes it in the work
        directory for the job
        '''
        slurm_params = self.calc_params.get('slurm_params', False)

        if slurm_params:
            txt = '#!/usr/bin/sh\n'
            for flag in  slurm_params.get('flags'):
                txt += f'#SBATCH {flag}\n'
            txt += '\n'
            for command in  slurm_params.get('precommands'):
                txt += command + '\n'
            txt += '\n'
            txt += 'echo "Job started on `hostname` at `date`"\n'
            txt += self._gen_run_command() + '\n'
            for command in slurm_params.get('postcommands'):
                txt += command + '\n'
            txt += 'echo "Job ended at `date`"'

        if write:
            slurm_file = join(self.workdir, 'job.sh')
            with open(slurm_file, 'w', encoding='utf-8') as f:
                f.write(txt)

    def _gen_slurm_interactive_txt(self):
        '''
        Generates an srun command to run the job
        '''
        slurm_params = self.calc_params.get('slurm_params', False)

        txt = 'srun '
        txt += ' '.join(slurm_params.get('flags'))
        txt += self._gen_run_command()

        return txt

    def gen_input_files(self):
        ''' Write input files to working directory '''
        workdir = self.workdir
        sim_params = self.sim_params.copy()
        calc_params = self.calc_params
        slurm = calc_params.get('slurm', False)
        batch = calc_params.get('batch', False)
        if os.path.exists(workdir) and not self.overwrite:
            print(f'Skipped writing in {self.workdir}')
        else:
            if not os.path.exists(workdir):
                os.makedirs(workdir)    
            cur_dir = os.getcwd()
            os.chdir(workdir)
            write_yaml('sim_params.yaml', sim_params)
            write_yaml('calc_params.yaml', calc_params)
            if slurm and batch:
                self._gen_slurm_script()

            os.chdir(cur_dir)

        self.ready = True


    def check_job_status(self):
        ''' Check job status'''
        if os.path.exists(self.output_file):
            output = read_yaml(self.output_file)
            status = output.get('status')
        else:
            status = 'not started'

        if status == 'complete':
            complete = True
        else:
            complete = False

        return complete, status

    def submit(self):
        ''' Submit a job using slurm, interactively or in the terminal '''

        if not self.ready:
            self.gen_input_files()

        if not self.sim_params.get('submit', False):
            return None

        cur_dir = os.getcwd()
        os.chdir(self.workdir)
        calc_params = self.calc_params
        sim_params = self.sim_params
        slurm = calc_params.get('slurm', False)
        interactive = calc_params.get('interactive', False)
        # Maybe want to do something clever with the exit code here
        if slurm and not interactive:
            command = 'sbatch job.sh'
        elif slurm and interactive:
            command = self._gen_slurm_interactive_txt().split()
        else:
            command = self._gen_run_command().split()

        try:
            completed_process = subprocess.run(
                command, check=True, capture_output=False
            )

            assert completed_process.returncode == 0, \
                f'Failed to run {sim_params.get("script")} in {self.workdir} \
                    with command {command}'
        except subprocess.CalledProcessError as exc:
            print("Status : FAIL", exc.returncode, exc.output)

        os.chdir(cur_dir)

# class ParallelJobs(Job):
#     ''' Unit job object with ability to submit independent jobs in parallel '''
#     # def __init__(
#     #     self,
#     #     sim_params: dict,
#     #     calc_params: dict,
#     #     workdir: str,
#     # ):
#     #     super(ParallelJobs, self).__init__(sim_params, calc_params, workdir)

#     def gen_slurm_array_script_txt(self, write: bool = True):
#         '''
#         Generates a slurm job array file and if necessary writes it in the work
#         directory for the job
#         '''
#         slurm_params = self.calc_params.get('slurm_params', False)
#         sim = self.sim_params.get('sim', False)

#         assert sim, 'Make sure you specify an existing script in sim_params'

#         if slurm_params:
#             txt = '#! /usr/bin/sh\n'
#             for command in  slurm_params.get('precommands'):
#                 txt += command + '\n'
#             txt += '\n'
#             txt += 'echo "Job started on `hostname` at `date`"\n'
#             txt += f'python {sim} sim_params.yaml calc_params.yaml\n'
#             for command in  slurm_params.get('postcommands'):
#                 txt += command + '\n'
#             txt += 'echo "Job Ended at `date`"'

#         slurm_script_txt = txt

#         if write:
#             slurm_file = join(workdir, 'job.sh')
#             with open(slurm_file, 'w', encoding='utf-8') as f:
#                 f.write(slurm_script_txt)

#     def submit_array(self):
#         if not self.slurm_script_txt:
#             self.gen_slurm_script_txt()
#         cur_dir = os.getcwd()
#         os.chdir(self.workdir)

#         # Maybe want to do something clever with the exit code here
#         subprocess.call('sbatch job.sh')
#         os.chdir(cur_dir)


def prepare_job(
    calc_params: dict,
    sim_params: dict,
    workdir: str,
):
    ''' Prepares a job's input files and gives a job object to handle '''
    job = UnitJob(sim_params, calc_params, workdir)
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
