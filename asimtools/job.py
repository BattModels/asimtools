'''
Set of tools for handling jobs in slurm, interactively or in the teminal for
a specific calculator.

Author: mkphuthi@github.com
'''

import subprocess
import os
from pathlib import Path
import functools
from glob import glob
from typing import List, TypeVar, Dict, Tuple, Union
from copy import deepcopy
import ase.io
import numpy as np
from asimtools.utils import (
    read_yaml,
    write_yaml,
    join_names,
    get_atoms,
)

Atoms = TypeVar('Atoms')

class Job():
    ''' Abstract class for the job object '''

    # pylint: disable=too-many-instance-attributes
    def __init__(
        self,
        config_input: Dict,
        sim_input: Dict,
    ):

        self.config_input = deepcopy(config_input)
        self.sim_input = deepcopy(sim_input)
        self.workdir = Path(
            sim_input.get('workdir', '.')
        ).resolve()
        self.launchdir = Path(os.getcwd())

        # Track if the job has associated image input file
        image = self.sim_input.get('args', {}).get('image', False)
        if image:
            self.atoms = get_atoms(**image)
        else:
            self.atoms = None

    def mkworkdir(self) -> None:
        ''' Creates the work directory if it doesn't exist '''
        if not self.workdir.exists():
            self.workdir.mkdir()

    def update_status(self, status: str) -> None:
        ''' Updates job status to specificied value '''
        self.update_output({'status': status})

    def start(self) -> None:
        ''' Updates the output to signal that the job was started '''
        self.update_status('started')

    def complete(self) -> None:
        ''' Updates the output to signal that the job was started '''
        self.update_status('complete')

    def fail(self) -> None:
        ''' Updates status to failed '''
        self.update_status('failed')

    def discard(self) -> None:
        ''' Updates status to discarded '''
        self.update_status('discard')

    def go_to_workdir(self) -> None:
        ''' Go to workdir '''
        if not self.workdir.exists():
            self.mkworkdir()
        os.chdir(self.workdir)

    def leave_workdir(self) -> None:
        ''' goes to directory from which job was launched '''
        os.chdir(self.launchdir)

    def add_output_files(self, file_dict: Dict) -> None:
        ''' Adds a file to the output file list '''
        files = self.get_output().get('files', {})
        files.update(file_dict)
        self.update_output({'files': file_dict})

    def set_input_image(self, atoms: Atoms) -> None:
        ''' Adds and writes the input image for simulation '''
        prefix = self.get_prefix()
        image_file = join_names([prefix, 'input_atoms.xyz'])
        self.sim_input['image'] = {
            'image_file': image_file
        }
        self.atoms = atoms

    def get_sim_input(self) -> Dict:
        ''' Get simulation input '''
        return deepcopy(self.sim_input)

    def get_prefix(self) -> str:
        ''' Get current job prefix '''
        return self.sim_input.get('prefix', '')

    def get_output_yaml(self) -> Dict:
        ''' Get current output file '''
        prefix = self.get_prefix()
        out_fname = join_names([prefix, 'output.yaml'])
        output_yaml = self.workdir / out_fname
        return output_yaml

    def get_config_input(self) -> Dict:
        ''' Get calculator input '''
        return deepcopy(self.config_input)

    def get_workdir(self) -> Path:
        ''' Get working directory '''
        return self.workdir

    def get_output(self) -> Dict:
        ''' Get values in output.yaml '''
        output_yaml = self.get_output_yaml()
        if output_yaml.exists():
            return read_yaml(output_yaml)
        else:
            return {}

    def get_output_images(self, ext='xyz') -> List:
        ''' Get all the output Atoms objects in workdir'''
        image_files = self.workdir / '*images_output.'+ext
        images = [ase.io.read(image_file) for image_file in image_files]
        assert len(images) > 0, f'No output images in {self.workdir}'

        return images

    def get_output_image(self, ext='extxyz') -> Atoms:
        ''' Get the output Atoms object '''
        files = self.get_output().get('files', {})
        image_file = files.get('image', False)
        if image_file:
            return ase.io.read(self.workdir / image_file, format=ext)
        else:
            return None

    def get_atoms(self) -> Atoms:
        ''' Get associated atoms object '''
        return self.atoms.copy()

    def update_sim_input(self, new_params) -> None:
        ''' Update simulation parameters '''
        self.sim_input.update(new_params)

    def update_calc_input(self, new_params) -> None:
        ''' Update calculator parameters '''
        self.config_input.update(new_params)

    def set_workdir(self, workdir) -> None:
        ''' Set working directory '''
        self.workdir = workdir

    def check_job_status(self) -> Tuple[bool,str]:
        ''' Check job status '''
        output = self.get_output()
        status = output.get('status', 'clean')
        complete = status == 'complete'
        return complete, status

    def update_output(self, output_update: Dict) -> None:
        ''' Update output.yaml if it exists or write a new one '''
        output = self.get_output()
        output.update(output_update)
        write_yaml(self.get_output_yaml(), output)

class UnitJob(Job):
    ''' Unit job object with ability to submit slurm jobs '''
    def __init__(
        self,
        config_input: Dict,
        sim_input: Dict,
    ):
        super().__init__(config_input, sim_input)

        self.script = sim_input.get('script', None)
        self.config_id = sim_input.get('config_id', None)
        self.calc_input = config_input[self.config_id]

    def gen_run_command(self) -> str:
        '''
        Generates the command to run the job, independent of slurm
        '''
        script = self.sim_input.get('script', False)
        assert script, 'Specify a script in sim_input'

        job_params = self.calc_input.get('job', {})
        run_prefix = job_params.get('run_prefix', '')
        run_suffix = job_params.get('run_suffix', '')

        txt = f'{run_prefix} '
        txt += 'asim-run calc_input.yaml sim_input.yaml '
        txt += f'{run_suffix} '
        return txt

    def _gen_slurm_script(self, write: bool = True) -> None:
        '''
        Generates a slurm job file and if necessary writes it in the work
        directory for the job
        '''
        slurm_params = self.calc_input.get('slurm', {})

        txt = '#!/usr/bin/sh\n'
        for flag in slurm_params.get('flags'):
            txt += f'#SBATCH {flag}\n'
        txt += '\n'
        for command in slurm_params.get('precommands', []):
            txt += command + '\n'
        txt += '\n'
        txt += 'echo "Job started on `hostname` at `date`"\n'
        txt += self.gen_run_command() + '\n'
        for command in slurm_params.get('postcommands', []):
            txt += command + '\n'
        txt += 'sleep 5\n'
        txt += 'echo "Job ended at `date`"'

        if write:
            slurm_file = self.workdir / 'job.sh'
            slurm_file.write_text(txt)

    def _gen_slurm_interactive_txt(self) -> str:
        '''
        Generates an srun command to run the job
        '''
        slurm_params = self.calc_input.get('slurm', False)

        txt = 'srun --unbuffered ' # To avoid IO lag
        txt += ' '.join(slurm_params.get('flags'))
        txt += self.gen_run_command()

        return txt

    def gen_input_files(self) -> None:
        ''' Write input files to working directory '''
        workdir = self.workdir
        sim_input = self.sim_input
        config_input = self.config_input
        calc_input = self.calc_input
        job_params = calc_input.get('job', {})
        use_slurm = job_params.get('use_slurm', False)
        interactive = job_params.get('interactive', False)
        overwrite = self.calc_input.get('overwrite', False)
        if workdir.exists() and not overwrite:
            print(f'Skipped writing in {self.workdir}')
        else:
            if not workdir.exists():
                workdir.mkdir()

            # # If we passed an atoms object, write to disk and link it
            # atoms = sim_input['args'].get('image', {}).get('atoms', False)
            # if atoms:
            #     image_file = str('image_input.xyz')
            #     atoms.write(self.get_workdir() / image_file)
            #     sim_input['args']['image']['image_file'] = image_file
            #     del sim_input['args']['image']['atoms']

            # Write the associated input image
            if self.atoms is not None:
                image_file = str('image_input.xyz')
                self.atoms.write(self.get_workdir() / image_file)
                sim_input['args']['image'] = {'image_file': image_file}

            images = self.sim_input['args'].get('images', False)
            if images:
                input_images_file = 'input_images.xyz'
                ase.io.write(input_images_file, images, format='extxyz')
                self.sim_input['args']['images'] = {
                    'image_file': input_images_file,
                }


            # # If we passed an atoms object, delete it otherwise
            # # the yamls can't store atoms objects
            # atoms = sim_input['args'].get('image', {}).get('atoms', False)
            # if atoms:
            #     del sim_input['args']['image']['atoms']

            write_yaml(workdir / 'sim_input.yaml', sim_input)
            # TODO: Consistency btn config and calc naming
            write_yaml(workdir / 'calc_input.yaml', config_input)
            write_yaml(self.get_output_yaml(), {'status': 'clean'})

        if use_slurm and not interactive:
            self._gen_slurm_script()

    def get_calc_input(self):
        ''' Get calc_input for this job specifically '''
        return deepcopy(self.calc_input)

    def submit(self, dependency: List = None) -> Union[None,List[str]]:
        ''' Submit a job using slurm, interactively or in the terminal '''

        self.gen_input_files()

        _, status = self.check_job_status()
        if status not in ('clean', 'discard'):
            return None

        if not self.calc_input.get('submit', False):
            return None

        cur_dir = Path('.').resolve()
        os.chdir(self.workdir)
        job_params = self.calc_input.get('job', {})
        sim_input = self.sim_input
        use_slurm = job_params.get('use_slurm', True)
        interactive = job_params.get('interactive', False)
        # Maybe want to do something clever with the exit code here
        if use_slurm and not interactive:
            if dependency is not None:
                dependstr = None
                for dep in dependency:
                    if dependstr is None:
                        dependstr = str(dep) if dep is not None else None
                    else:
                        dependstr += f':{dep}' if dep is not None else ''

                command = ['sbatch', '-d', f'afterok:{dependstr}', 'job.sh']
            else:
                command = ['sbatch', 'job.sh']
        elif use_slurm and interactive:
            command = self._gen_slurm_interactive_txt().split()
        else:
            command = self.gen_run_command().split()

        run_job = False
        overwrite = sim_input.get('overwrite', False)
        _, status = self.check_job_status()
        if overwrite or status == 'clean':
            run_job = True

        if run_job:
            completed_process = subprocess.run(
                command, check=False, capture_output=True, text=True,
            )

            with open('stdout.txt', 'w', encoding='utf-8') as output_file:
                output_file.write(completed_process.stdout)

            if completed_process.returncode != 0:
                err_txt = f'Failed to run {sim_input.get("script")} in '
                err_txt += f'{self.workdir} with command:\n'
                err_txt += f'{" ".join(command)}\n'
                err_txt += f'See {self.workdir / "stderr.txt"} for details.'
                print(err_txt)
                with open('stderr.txt', 'w', encoding='utf-8') as err_file:
                    err_file.write(completed_process.stderr)
                completed_process.check_returncode()

        os.chdir(cur_dir)

        if use_slurm:
            job_ids = [int(completed_process.stdout.split(' ')[-1])]
            return job_ids
        return None

class DistributedJob(Job):
    ''' Array job object with ability to submit simultaneous jobs '''
    def __init__(
        self,
        config_input: Dict,
        sim_input: Dict,
    ):
        super().__init__(config_input, sim_input)

        # self.script = sim_input.get('script', None)
        # self.config_id = sim_input.get('config_id', None)
        # self.calc_input = config_input[self.config_id]

        # Set a standard for all subdirectories to start
        # with "id-{job_id}" followed by user labels
        new_sim_input = {}
        for i, (sim_id, subsim_input) in enumerate(sim_input.items()):
            assert 'script' in subsim_input, 'Check sim_input format,\
                must have script and each job with a unique key'

            if not sim_id.startswith('id-'):
                new_sim_id = join_names([f'id-{i}',sim_id])
                new_sim_input[new_sim_id] = deepcopy(sim_input[sim_id])

        sim_input = new_sim_input
        unitjobs = []
        for sim_id, subsim_input in sim_input.items():
            assert 'script' in subsim_input, 'Check sim_input format,\
                must have script and each job with a unique key'

            if not subsim_input.get('workdir', False):
                subsim_input['workdir'] = sim_id
            unitjob = UnitJob(config_input, subsim_input)
            unitjobs.append(unitjob)

        # If all the jobs have the same config, use a job array
        config_id = unitjobs[0].config_id
        config0 = config_input[config_id]
        same_config = np.all(
            [(uj.config_id == config_id) for uj in unitjobs]
        )

        if same_config and config0['job'].get('use_slurm', False):
            self.use_array = True
        else:
            self.use_array = False

        self.unitjobs = unitjobs

    def submit_jobs(self) -> Union[None,List[str]]:
        '''
        Submits the jobs. If submitting lots of batch jobs, we 
        recommend using DistributedJob.submit_array
        '''
        job_ids = []
        for unitjob in self.unitjobs:
            job_id = unitjob.submit()
            job_ids.append(job_id)
        return job_ids

    def submit_array(
        self,
        array_max=None,
        dependency: List = None
    ) -> Union[None,List[str]]:
        '''
        Submits a job array if all the jobs have the same config_id
        '''
        self.gen_input_files()

        njobs = len(self.unitjobs)

        arr_max_str = ''
        if array_max is not None:
            arr_max_str = f'%{array_max}'

        # command = f'sbatch --array=[0-{njobs-1}] job_array.sh'
        if dependency is not None:
            dependstr = None
            for dep in dependency:
                if dependstr is None:
                    dependstr = str(dep) if dep is not None else None
                else:
                    dependstr += f':{dep}' if dep is not None else ''

            command = [
                'sbatch',
                f'--array=[0-{njobs-1}]{arr_max_str}',
                '-d', f'afterok:{dependstr}',
                'job_array.sh'
            ]
        else:
            command = [
                'sbatch',
                f'--array=[0-{njobs-1}]{arr_max_str}',
                'job_array.sh'
            ]

        completed_process = subprocess.run(
            command, check=False, capture_output=True, text=True,
        )

        with open('stdout.txt', 'w', encoding='utf-8') as output_file:
            output_file.write(completed_process.stdout)

        if completed_process.returncode != 0:
            err_txt = 'Failed to run array in '
            err_txt += f'{self.workdir} with command:\n'
            err_txt += f'{" ".join(command)}\n'
            err_txt += f'See {self.workdir / "stderr.txt"} for details.'
            print(err_txt)
            with open('stderr.txt', 'w', encoding='utf-8') as err_file:
                err_file.write(completed_process.stderr)
            completed_process.check_returncode()

        job_ids = [int(completed_process.stdout.split(' ')[-1])]
        return job_ids

    def _gen_array_script(self, write: bool = True) -> None:
        '''
        Generates a slurm job array file and if necessary 
        writes it in the work directory for the job. Only works
        if there is one config_id for all jobs
        '''
        config_id = self.unitjobs[0].config_id
        config = self.config_input[config_id]

        slurm_params = config.get('slurm', {})

        txt = '#!/usr/bin/sh\n'
        for flag in slurm_params.get('flags'):
            txt += f'#SBATCH {flag}\n'

        txt += 'if [[ ! -z ${SLURM_ARRAY_TASK_ID} ]]; then\n'
        txt += '    fls=( id-* )\n'
        txt += '    WORKDIR=${fls[${SLURM_ARRAY_TASK_ID}]}\n'
        txt += 'fi\n\n'
        txt += 'CUR_DIR=`pwd`\n'
        txt += 'cd ${WORKDIR}\n'
        txt += '\n'
        txt += '\n'.join(slurm_params.get('precommands', []))
        txt += '\n'
        txt += 'echo "LAUNCHDIR: and ${CUR_DIR}"\n'
        txt += 'echo "WORKDIR: ${WORKDIR}"\n'
        txt += 'echo "Job started on `hostname` at `date`"\n'
        txt += self.unitjobs[0].gen_run_command() + '\n'
        txt += '\n'.join(slurm_params.get('postcommands', []))
        txt += '\n'
        txt += 'cd ${CUR_DIR}\n'
        # Wait a few seconds for the so that files written to disk
        # don't lag. Might be worth switching to sbatch --unbuffered
        txt += 'sleep 3\n'
        txt += 'echo "Job ended at `date`"'

        if write:
            slurm_file = self.workdir / 'job_array.sh'
            slurm_file.write_text(txt)

    def gen_input_files(self) -> None:
        ''' Write input files to working directory '''
        # # If we passed an atoms object or images, delete. Otherwise
        # # the yamls can't store atoms objects
        # images = self.sim_input['args'].get('images', False)
        # if images:
        #     ase.io.write('input_images.xyz', images, format='extxyz')
        #     del self.sim_input['args']['images']['images']

        if self.use_array:
            self._gen_array_script()
        for unitjob in self.unitjobs:
            unitjob.gen_input_files()

    def submit(self) -> None:
        ''' 
        Submit a job using slurm, interactively or in the terminal
        '''

        if not self.workdir.exists():
            self.gen_input_files()

        _, status = self.check_job_status()

        if status in ('completed'):
            if not self.sim_input.get('overwrite', False):
                return None

        if not self.sim_input.get('submit', True):
            return None
        cur_dir = Path('.').resolve()
        os.chdir(self.workdir)

        job_ids = None
        if self.use_array:
            job_ids = self.submit_array()
        else:
            job_ids = self.submit_jobs()

        os.chdir(cur_dir)
        return job_ids

class ChainedJob(Job):
    ''' Jobs made of smaller sequential jobs '''
    def __init__(
        self,
        config_input: Dict,
        sim_input: Dict,
    ):
        super().__init__(config_input, sim_input)

        assert np.all([key.startswith('step-') for key in sim_input]),\
            'Keys should take the form "step-{step_id}" from 1'

        # workdir naming should follow pattern. Might loosen this
        # eventually
        for sim_id, subsim_input in sim_input.items():
            subsim_input['workdir'] = sim_id

        unitjobs = []
        for i in range(len(sim_input)):
            subsim_input = deepcopy(sim_input[f'step-{i}'])
            unitjob = UnitJob(config_input, subsim_input)
            unitjobs.append(unitjob)

        # # If all the jobs have the same config, use a job array
        # config_id = unitjobs[0].config_id
        # config0 = config_input[config_id]
        # same_config = np.all(
        #     [(uj.config_id == config_id) for uj in unitjobs]
        # )

        # if same_config and config0['job'].get('use_slurm', False):
        #     self.use_array = True
        # else:
        #     self.use_array = False

        self.unitjobs = unitjobs

    def get_last_output(self) -> Dict:
        ''' Returns the output of the last job in the chain '''
        return self.unitjobs[-1].get_output()

    def get_current_step(self):
        ''' Gets current step in the chain '''
        raise NotImplementedError

    def submit(self, dependency: List = None) -> List:
        ''' 
        Submit a job using slurm, interactively or in the terminal
        '''

        # if not self.workdir.exists():
        #     self.mkworkdir()

        # _, status = self.check_job_status()

        # if status in ('completed'):
        #     if not self.sim_input.get('overwrite', False):
        #         return None

        # if not self.sim_input.get('submit', True):
        #     return None

        cur_dir = Path('.').resolve()
        os.chdir(self.workdir)

        step = 0 #self.get_current_step()
        are_interactive_jobs = [
            uj.get_calc_input()['job'].get('interactive', False) \
                for uj in self.unitjobs
        ]
        use_batch = not np.any(are_interactive_jobs)
        if use_batch:
            if step != len(self.unitjobs):
                dependency = None
                for unitjob in self.unitjobs[step:]:
                    dependency = unitjob.submit(dependency=dependency)

            job_ids = dependency
        else:
            for unitjob in self.unitjobs[step:]:
                job_ids = unitjob.submit()

        os.chdir(cur_dir)

        return job_ids

# def check_jobs(jobs) -> bool:
#     ''' Checks status of jobs and prints a report '''

#     completed = []
#     print('+'*20)
#     print('Job statuses: ')
#     print('+'*20)
#     for job in jobs:
#         complete, status = job.check_job_status()
#         print(f'{status}: {job.get_workdir()}')
#         if complete or status == 'discard':
#             completed.append(True)
#         else:
#             completed.append(False)

#     if np.all(completed):
#         return True

#     return False

def leaf(func):
    ''' 
    Step wrapper. This handles things like:
    - Going to the correct workdir
    - initializing output files for progress tracking
    while being minimally invasive
    '''
    # functools makes it so that the docs work instead of giving the
    # wrapper docstring
    @functools.wraps(func)
    def wrapper(*arg, **kwargs):
        workdir = kwargs.get('workdir', '.')
        prefix = kwargs.get('prefix', '')
        config_id = kwargs.get('config_id', False)
        assert config_id, 'No config_id in leaf'
        sim_input = {
            'workdir': workdir,
            'prefix': prefix,
        }
        config_input = arg[0]
        calc_input = config_input[config_id]
        job = Job(*arg, sim_input=sim_input)
        job.go_to_workdir()
        job.start()
        try:
            results = func(calc_input, *arg[1:], **kwargs)
        except:
            job.fail()
            raise
        job.update_output(results)
        job.complete()
        job.leave_workdir()
        return results
    return wrapper

def branch(func):
    ''' 
    Step wrapper. This handles things like:
    - Going to the correct workdir
    - initializing output files for progress tracking
    while being minimally invasive
    '''
    # functools makes it so that the docs work instead of giving the
    # wrapper docstring
    @functools.wraps(func)
    def wrapper(*arg, **kwargs):
        workdir = kwargs.get('workdir', '.')
        prefix = kwargs.get('prefix', '')
        sim_input = {
            'workdir': workdir,
            'prefix': prefix,
        }
        job = Job(*arg, sim_input=sim_input)
        job.go_to_workdir()
        job.start()
        try:
            results = func(*arg, **kwargs)
        except:
            job.fail()
            raise
        job.update_output(results)
        job.complete()
        job.leave_workdir()
        return results
    return wrapper


def load_jobs_from_directory(
    workdir: str = None,
    workdirs: list = None,
    pattern: str = None,
) -> List[UnitJob]:
    ''' Loads all jobs from a given directory/directories '''

    assert [workdir, workdirs, pattern].count(None) == 2, \
        'Provide exactly one of workdir, workdirs or pattern'

    job_list = []
    if pattern is not None:
        job_list += load_jobs_from_directory(glob(pattern))
    elif workdirs is not None:
        for wd in workdirs:
            job_list += load_jobs_from_directory(wd)
    else:
        calc_input = read_yaml(
            os.path.join(workdir, 'calc_input.yaml')
        )
        sim_input = read_yaml(
            os.path.join(workdir, 'sim_input.yaml')
        )
        job_list += [UnitJob(calc_input, sim_input)]

    return job_list

def load_outputs(
    workdir: str = None,
    workdirs: list = None,
    pattern: str = None,
) -> List[UnitJob]:
    ''' Loads all jobs from a given directory/directories '''

    assert [workdir, workdirs, pattern].count(None) == 2, \
        'Provide exactly one of workdir, workdirs or pattern'

    if workdirs is not None or pattern is not None:
        outputs = []
        if pattern is not None:
            workdirs = glob(pattern)

        for workdir in workdirs:
            outputs += load_outputs(workdir=workdir)
    else:
        output_files = glob(str(Path(workdir) / '*output.yaml'))
        assert len(output_files) == 1, \
            f'{len(output_files)} output files in {workdir}'

        outputs = [{
            'path': str(workdir),
            'output': read_yaml(output_files[0])
        }]

    return outputs

def load_output_property(
    prop: str,
    workdir: str = None,
    workdirs: list = None,
    pattern: str = None,
) -> List[UnitJob]:
    ''' Loads all jobs from a given directory/directories '''

    assert [workdir, workdirs, pattern].count(None) == 2, \
        'Provide exactly one of workdir, workdirs or pattern'

    outputs = load_outputs(
        workdir=workdir,
        workdirs=workdirs,
        pattern=pattern
    )
    properties = []
    for output in outputs:
        properties.append(output[prop])

    return properties

# def load_output_images(
#     workdir: str = None,
#     workdirs: list = None,
#     pattern: str = None,
# ) -> List[UnitJob]:
#     ''' Loads all jobs from a given directory/directories '''

#     output_files = load_output_property(
#         'files',
#         workdir=workdir,
#         workdirs=workdirs,
#         pattern=pattern
#     )

#     image_files = [outf['image'] for outf in output_files]
#     output_images = []
#     for image_file in image_files:
#         output_images.append(ase.io.read(image_file))

#     return output_images

def load_input_images(
    workdir: str = None,
    workdirs: list = None,
    pattern: str = None,
) -> List[UnitJob]:
    ''' Loads all jobs from a given directory/directories '''

    assert [workdir, workdirs, pattern].count(None) == 2, \
        'Provide exactly one of workdir, workdirs or pattern'

    if workdirs is not None or pattern is not None:
        if pattern is not None:
            workdirs = glob(pattern)

        for workdir in workdirs:
            input_images += load_input_images(workdir=workdir)
    else:
        input_images = glob(str(Path(workdir) / 'image_input.xyz'))
        assert len(input_images) == 1, \
            f'{len(input_images)} input image files in {workdir}'

        input_images = [ase.io.read(input_images[0])]

    return input_images

def load_output_images(
    workdir: str = None,
    workdirs: list = None,
    pattern: str = None,
) -> List[UnitJob]:
    ''' Loads all jobs from a given directory/directories '''

    assert [workdir, workdirs, pattern].count(None) == 2, \
        'Provide exactly one of workdir, workdirs or pattern'

    if workdirs is not None or pattern is not None:
        if pattern is not None:
            workdirs = glob(pattern)

        output_images = []
        for workdir in workdirs:
            output_images += load_output_images(workdir=workdir)
    else:
        output_images = glob(str(Path(workdir) / 'image_output.xyz'))
        import os; print('CWD', os.getcwd())
        assert len(output_images) == 1, \
            f'{len(output_images)} output image files in {workdir}'

        output_images = [ase.io.read(output_images[0])]

    return output_images

# def prepare_job(
#     calc_input: Dict,
#     sim_input: Dict,
#     workdir: str,
# ) -> Job:
#     ''' Prepares a job's input files and gives a job object to handle '''
#     job = UnitJob(calc_input, sim_input, workdir)
#     job.gen_input_files()
#     return job


# def gen_slurm_array_txt(calc_input: Dict, sim_input: Dict) -> str:
#     '''
#     Generates the job array script from given data following standard pattern 
#     '''
#     slurm_params = calc_input['slurm_params']
#     script = sim_input['script']
#     flags = slurm_params['flags']
#     precommands = slurm_params['precommands']
#     postcommands = slurm_params['postcommands']
#     txt = '#!/usr/bin/bash\n'
#     txt += '\n'.join([f'#SBATCH {flag}' for flag in flags])
#     txt += '\n'
#     txt += 'if [[ ! -z ${SLURM_ARRAY_TASK_ID} ]]; then\n'
#     txt += '    fls=( id_* )\n'
#     txt += '    WORKDIR=${fls[${SLURM_ARRAY_TASK_ID}]}\n'
#     txt += 'fi\n\n'
#     txt += 'cd ${WORKDIR}\n'
#     txt += 'echo "Job started on `hostname` at `date`"\n'
#     txt += '\n'.join(precommands)
#     txt += f'{script} calc_input sim_input'
#     txt += '\n'.join(postcommands)
#     txt += 'echo "Job Ended at `date`"'

#     return txt

# def prepare_array_jobs(jobs):
#     ''' Prepares multiple jobs and a job array for them '''
#     workdirs = []
#     prefixes = []
#     for job in jobs:
#         job.gen_input_files()
#         workdirs.append(job.get_workdir())

#     slurm_params = job.get_calc_input()['slurm']
#     script = job.get_sim_input()['script']
#     gen_slurm_array_txt(script, slurm_params, workdirs, prefixes)
