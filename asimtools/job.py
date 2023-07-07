'''
Set of tools for handling jobs in slurm, interactively or in the teminal for
a specific calculator.

Author: mkphuthi@github.com
'''

import subprocess
import os
from pathlib import Path
from datetime import datetime
# import functools
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
    get_images,
    get_env_input,
    get_calc_input,
    check_if_slurm_job_is_running,
)

Atoms = TypeVar('Atoms')

class Job():
    ''' Abstract class for the job object '''

    # pylint: disable=too-many-instance-attributes
    def __init__(
        self,
        sim_input: Dict,
        env_input: Union[Dict,None] = None,
        calc_input: Union[Dict,None] = None,
    ) -> None:
        if env_input is None:
            env_input = get_env_input()
        if calc_input is None:
            calc_input = get_calc_input()
        self.env_input = deepcopy(env_input)
        self.calc_input = deepcopy(calc_input)
        self.sim_input = deepcopy(sim_input)
        self.workdir = Path(
            sim_input.get('workdir', '.')
        ).resolve()
        # Any subsequent input files should use the full path
        self.sim_input['workdir'] = str(self.workdir)
        self.launchdir = Path(os.getcwd())

    def mkworkdir(self) -> None:
        ''' Creates the work directory if it doesn't exist '''
        if not self.workdir.exists():
            self.workdir.mkdir()

    def update_status(self, status: str) -> None:
        ''' Updates job status to specificied value '''
        self.update_output({'status': status})

    def start(self) -> None:
        ''' Updates the output to signal that the job was started '''
        self.update_output({'start_time': datetime.now()})
        self.update_status('started')
        print(self.get_output())

    def complete(self) -> None:
        ''' Updates the output to signal that the job was started '''
        self.update_output({'end_time': datetime.now()})
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

    def get_sim_input(self) -> Dict:
        ''' Get simulation input '''
        return deepcopy(self.sim_input)

    def get_calc_input(self) -> Dict:
        ''' Get calculator input '''
        return deepcopy(self.calc_input)

    def get_env_input(self) -> Dict:
        ''' Get environment input '''
        return deepcopy(self.env_input)

    # def get_prefix(self) -> str:
    #     ''' Get current job prefix '''
    #     return self.sim_input.get('prefix', '')

    def get_output_yaml(self) -> Dict:
        ''' Get current output file '''
        # prefix = self.get_prefix()
        # Will experiment with prefixes once everything works
        prefix = ''
        out_fname = join_names([prefix, 'output.yaml'])
        output_yaml = self.workdir / out_fname
        return output_yaml

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

    # def get_output_images(self, ext='xyz') -> List:
    #     ''' Get all the output Atoms objects in workdir'''
    #     files = self.get_output().get('files', {})
    #     image_file = files.get('images', False)
    #     if image_file:
    #         return ase.io.read(
    #             self.workdir / image_file, format=ext, index=':',
    #         )
    #     else:
    #         return None

    # def get_output_image(self, ext='extxyz') -> Union[Atoms,None]:
    #     ''' Get the output Atoms object '''
    #     files = self.get_output().get('files', {})
    #     image_file = files.get('image', False)
    #     if image_file:
    #         return ase.io.read(self.workdir / image_file, format=ext)
    #     else:
    #         return None

    def update_sim_input(self, new_params) -> None:
        ''' Update simulation parameters '''
        self.sim_input.update(new_params)

    def update_calc_input(self, new_params) -> None:
        ''' Update calculator parameters '''
        self.calc_input.update(new_params)

    def update_env_input(self, new_params) -> None:
        ''' Update calculator parameters '''
        self.env_input.update(new_params)

    def set_workdir(self, workdir) -> None:
        ''' Set working directory '''
        self.workdir = workdir

    # def set_atoms(self, atoms) -> None:
    #     ''' Set atoms object '''
    #     self.atoms = atoms

    def get_status(self, display=False) -> Tuple[bool,str]:
        ''' Check job status '''
        output = self.get_output()
        job_id = output.get('job_id', False)
        running = False
        if job_id:
            running = check_if_slurm_job_is_running(job_id)
        if running:
            status = 'started'
            self.update_status(status)
            complete = False
        else:
            status = output.get('status', 'clean')
            complete = status == 'complete'
        if display:
            print(f'{status}: {self.workdir}')
        return complete, status

    def update_output(self, output_update: Dict) -> None:
        ''' Update output.yaml if it exists or write a new one '''
        output = self.get_output()
        output.update(output_update)
        write_yaml(self.get_output_yaml(), output)

class UnitJob(Job):
    ''' 
    Unit job object with ability to submit a slurm/interactive job.
    Unit jobs run in a specific environment and if required, run a 
    specific calculator. More complex workflows are built up of unit
    jobs
    '''
    def __init__(
        self,
        sim_input: Dict,
        env_input: Union[Dict,None] = None,
        calc_input: Union[Dict,None] = None,
    ) -> None:
        super().__init__(sim_input, env_input, calc_input)
        self.env_id = self.sim_input.get('env_id', None)
        if self.env_id is not None:
            self.env = self.env_input[self.env_id]
        else:
            self.env = {'mode': {
                    'use_slurm': False,
                    'interactive': True,
                }
            }

        # print('unitjob init', calc_input)
        self.calc_id = self.sim_input.get('calc_id', None)
        if self.calc_id is not None:
            self.calc_params = self.calc_input[self.calc_id]

    def gen_run_command(self) -> str:
        '''
        Generates the command to run the job, independent of slurm
        '''
        script = self.sim_input.get('script', False)
        assert script, 'Specify a script in sim_input'

        mode_params = self.env.get('mode', {
            'use_slurm': False, 'interactive': True
        })
        run_prefix = mode_params.get('run_prefix', '')
        run_suffix = mode_params.get('run_suffix', '')

        txt = f'{run_prefix} '
        txt += 'asim-run sim_input.yaml -c calc_input.yaml'
        txt += f'{run_suffix} '
        return txt

    def _gen_slurm_script(self, write: bool = True) -> None:
        '''
        Generates a slurm job file and if necessary writes it in 
        the work directory for the job
        '''
        slurm_params = self.env.get('slurm', {})

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
        # txt += 'sleep 5\n'
        txt += 'echo "Job ended at `date`"'

        if write:
            slurm_file = self.workdir / 'job.sh'
            slurm_file.write_text(txt)

    def _gen_slurm_interactive_txt(self) -> str:
        '''
        Generates an srun command to run the job
        '''
        slurm_params = self.env.get('slurm', False)

        txt = 'salloc --unbuffered ' # unbuffered to avoid IO lag
        txt += ' '.join(slurm_params.get('flags'))
        txt += self.gen_run_command()

        return txt

    def gen_input_files(self) -> None:
        ''' Write input files to working directory '''
        workdir = self.workdir
        sim_input = deepcopy(self.sim_input)
        # The sim_input file will already be in the work directory
        sim_input['workdir'] = '.'
        mode_params = self.env.get('mode', {})
        use_slurm = mode_params.get('use_slurm', False)
        interactive = mode_params.get('interactive', False)
        overwrite = self.sim_input.get('overwrite', False)
        output_file = workdir / 'output.yaml'
        sim_input_file = workdir / 'sim_input.yaml'
        if output_file.exists() and sim_input_file.exists():
            started = read_yaml(output_file).get('start_time', False)
            if started and not overwrite:
                txt = f'Skipped writing in {self.workdir}, job already '
                txt += 'started. Set overwrite=True to overwrite'
                print(txt)
                return None
        if not workdir.exists():
            workdir.mkdir()

        image = self.sim_input.get('args', {}).get('image', False)

        if image:
            # if image.get('atoms', False):
            atoms = get_atoms(**image)
            input_image_file = self.workdir / 'input_image.xyz'
            atoms.write(input_image_file, format='extxyz')
            self.sim_input['args']['image'] = {
                'image_file': str(input_image_file),
            }

        images = self.sim_input.get('args', {}).get('images', False)
        if images:
            if images.get('images', False):
                images = get_images(**images)
                input_images_file = self.workdir / 'input_images.xyz'
                # print('writing images', image, images)
                ase.io.write(input_images_file, images, format='extxyz')
                self.sim_input['args']['images'] = {
                    'image_file': str(input_images_file),
                }
                # print('new image', self.sim_input)
            elif images.get('image_file', False):
                # import os
                # print('image_file cwd', os.getcwd())
                # print('image_file before', images['image_file'])
                input_images_file = Path(images['image_file']).resolve()
                self.sim_input['args']['images'] = {
                    'image_file': str(input_images_file),
                }
                # print('image_file after', input_images_file)

        write_yaml(workdir / 'sim_input.yaml', self.sim_input)
        write_yaml(workdir / 'calc_input.yaml', self.calc_input)
        write_yaml(workdir / 'env_input.yaml', self.env_input)
        write_yaml(self.get_output_yaml(), {'status': 'clean'})

        if use_slurm and not interactive:
            self._gen_slurm_script()
        return None

    def submit(
        self, dependency: Union[List,None] = None
    ) -> Union[None,List[str]]:
        '''
        Submit a job using slurm, interactively or in the terminal
        '''
        self.gen_input_files()
        # print('unitjob submit: start')
        _, status = self.get_status(display=True)
        if status in ('started', 'discard'):
            return None

        if not self.sim_input.get('submit', False):
            return None

        print('unitjob submit: overwritten files')
        print(f'unitjob: submit: {self.sim_input.get("submit")}')
        cur_dir = Path('.').resolve()
        os.chdir(self.workdir)
        mode_params = self.env.get('mode', {})
        # print('unitjob mode_params:', mode_params)
        sim_input = self.sim_input
        use_slurm = mode_params.get('use_slurm', False)
        interactive = mode_params.get('interactive', True)
        # print('unitjob submit: setup parameters')
        if use_slurm and not interactive:
            # print('unitjob submit: using slurm')
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
            # print('unitjob submit: interactive slurm')
            command = self._gen_slurm_interactive_txt().split()
        else:
            # print('unitjob submit: inline')
            command = self.gen_run_command().split()

        run_job = False
        overwrite = sim_input.get('overwrite', False)
        _, status = self.get_status()
        if overwrite or status == 'clean':
            run_job = True

        # print(f'unitjob submit: overwrite: {overwrite}, status: {status}')
        # print('unitjob submit: runjob: ', run_job)
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
        # print('unitjob submit: finished job')
        if use_slurm and not interactive:
            # print('completed_process:', completed_process.stdout)
            job_ids = [int(completed_process.stdout.split(' ')[-1])]
            # print('unitjob submit: job_ids:', job_ids)
            return job_ids

        return None

class DistributedJob(Job):
    ''' Array job object with ability to submit simultaneous jobs '''
    def __init__(
        self,
        sim_input: Dict,
        env_input: Union[None,Dict] = None,
        calc_input: Union[None,Dict] = None,
    ) -> None:
        super().__init__(sim_input, env_input, calc_input)
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

            # if not subsim_input.get('workdir', False):
            subsim_input['workdir'] = sim_id
            unitjob = UnitJob(subsim_input, env_input, calc_input)
            unitjobs.append(unitjob)

        # If all the jobs have the same config, use a job array
        env_id = unitjobs[0].env_id
        # env0 = env_input[env_id]
        same_env = np.all(
            [(uj.env_id == env_id) for uj in unitjobs]
        )
        # calc_id = unitjobs[0].calc_id
        # calc0 = calc_input[calc_id]
        # same_calc = np.all(
        #     [(uj.calc_id == calc_id) for uj in unitjobs]
        # )
        all_slurm = np.all(
            [uj.env['mode'].get('use_slurm', False) for uj in unitjobs]
        )

        if same_env and all_slurm:
            self.use_array = True
        else:
            self.use_array = False

        self.unitjobs = unitjobs

    def submit_jobs(self, **kwargs) -> Union[None,List[str]]:
        '''
        Submits the jobs. If submitting lots of batch jobs, we 
        recommend using DistributedJob.submit_array
        '''
        # It is possible to do an implementation where we limit
        # the number of simultaneous jobs like an array using job
        # dependencies but that can be implemented later
        job_ids = []
        for unitjob in self.unitjobs:
            job_id = unitjob.submit()
            job_ids.append(job_id)
        return job_ids

    def submit_array(
        self,
        array_max=None,
        dependency: Union[List[str],None] = None,
        **kwargs,
    ) -> Union[None,List[str]]:
        '''
        Submits a job array if all the jobs have the same env and use slurm
        '''
        self.gen_input_files()

        njobs = len(self.unitjobs)

        # For limiting number of jobs launched in array
        if array_max is not None:
            arr_max_str = f'%{array_max}'
        else:
            arr_max = self.unitjobs[0].env['mode'].get('array_max', False)
            if arr_max:
                arr_max_str = f'%{arr_max}'
            else:
                arr_max_str = ''


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
        env_id = self.unitjobs[0].env_id
        env = self.env_input[env_id]

        slurm_params = env.get('slurm', {})

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
        # txt += 'sleep 3\n'
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

    def submit(self, **kwargs) -> None:
        ''' 
        Submit a job using slurm, interactively or in the terminal
        '''
        # print('djob: submit')
        # if not self.workdir.exists():
        #     self.gen_input_files()
        # self.gen_input_files()
        # _, status = self.check_job_status()

        # if status in ('completed'):
        #     if not self.sim_input.get('overwrite', False):
        #         return None

        # if not self.sim_input.get('submit', True):
        #     return None
        cur_dir = Path('.').resolve()
        os.chdir(self.workdir)
        job_ids = None
        if self.use_array:
            # print('djob: submit array')
            job_ids = self.submit_array(**kwargs)
        else:
            # print('djob: submit individual jobs')
            job_ids = self.submit_jobs(**kwargs)
        # print('djob job_ids', job_ids)
        os.chdir(cur_dir)
        return job_ids

class ChainedJob(Job):
    '''
    Jobs made of smaller sequential unitjobs. Note that this only works
    well with unit jobs. If one of the UnitJobs calls scripts, internally,
    slurm will fail to build the correct dependencies but the files will be
    generated appropriately, you can submit them manually
    '''
    def __init__(
        self,
        sim_input: Dict,
        env_input: Union[None,Dict] = None,
        calc_input: Union[None,Dict] = None,
    ) -> None:
        super().__init__(sim_input, env_input, calc_input)

        assert np.all([key.startswith('step-') for key in sim_input]),\
            'Keys should take the form "step-{step_id}" from step_id=0'

        # workdir naming should follow pattern. Might loosen this
        # eventually
        for sim_id, subsim_input in sim_input.items():
            subsim_input['workdir'] = sim_id

        unitjobs = []
        for i in range(len(sim_input)):
            subsim_input = deepcopy(sim_input[f'step-{i}'])
            unitjob = UnitJob(subsim_input, env_input, calc_input)
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

    def get_current_step(self) -> int:
        ''' Gets current step in the chain '''
        raise NotImplementedError

    def submit(self, dependency: Union[List,None] = None) -> List:
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
            uj.env['mode'].get('interactive', False) \
                for uj in self.unitjobs
        ]
        all_interactive = np.all(are_interactive_jobs)
        all_not_interactive = np.sum(are_interactive_jobs) == 0
        assert all_interactive or all_not_interactive, \
            'All jobs in chain must be interactive or not interactive'

        are_using_slurm = [
            uj.env['mode'].get('use_slurm', False) \
                for uj in self.unitjobs
        ]
        all_slurm = np.all(are_using_slurm)
        all_not_slurm = np.sum(are_using_slurm) == 0

        assert all_slurm or all_not_slurm, \
            'All jobs in chain must use slurm or none'

        # Only use slurm depndencies if all the jobs are batch jobs
        if all_slurm and all_not_interactive:
            if step != len(self.unitjobs):
                dependency = None
                for i, unitjob in enumerate(self.unitjobs[step:]):
                    # if step + i > 0:
                    #     prev_step_complete, _ = self.unitjobs[step:][i-1].get_status()
                    # else:
                    #     prev_step_complete = True
                    # print(f'step: {step+i}, prev_complete: {prev_step_complete}')
                    # if prev_step_complete:
                    cur_step_complete, status = unitjob.get_status()
                    print(f'step: {step+i}, cur_complete: {cur_step_complete}')
                    if not cur_step_complete and status != 'started':
                        unitjob.env['slurm']['flags'].append(f'-J step-{step+i}')
                        # print('Dependency:', dependency, unitjob.sim_input)
                        dependency = unitjob.submit(dependency=dependency)
                        # print('Dependency uj output:', unitjob.get_output())
                        # dependency = unitjob.get_output().get('job_ids', None)
                        # print('Dependency next:', dependency)

            job_ids = dependency
        else:
            for unitjob in self.unitjobs[step:]:
                job_ids = unitjob.submit()

        os.chdir(cur_dir)

        return job_ids


def load_job_from_directory(workdir: str):
    ''' Loads a job from a given directory '''
    workdir = Path(workdir)
    sim_inputs = glob(str(workdir / 'sim_input.yaml'))
    if len(sim_inputs) != 1:
        print(f'WARNING: Multiple or no sim_input files in {workdir}')
    try:
        sim_input = read_yaml(glob(str(workdir / 'sim_input.yaml'))[0])
    except:
        print('sim_input.yaml not found')
        raise

    env_inputs = glob(str(workdir / 'env_input.yaml'))
    if len(env_inputs) == 1:
        env_input = read_yaml(env_inputs[0])
    else:
        env_input = None

    calc_inputs = glob(str(workdir / 'calc_input.yaml'))
    if len(calc_inputs) == 1:
        calc_input = read_yaml(calc_inputs[0])
    else:
        calc_input = None

    job = Job(sim_input, env_input, calc_input)
    return job


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
        assert len(output_images) == 1, \
            f'{len(output_images)} output image files in {workdir}'

        output_images = [ase.io.read(output_images[0])]

    return output_images

################ Graveyard of failed attempts ##########################
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

# def wrap_job(func):
#     ''' 
#     Step wrapper. This handles things like:
#     - Going to the correct workdir
#     - initializing output files for progress tracking
#     while being minimally invasive
#     '''
#     # functools makes it so that the docs work instead of giving the
#     # wrapper docstring
#     @functools.wraps(func)
#     def wrapper(*args, **kwargs):
#         sim_input = {}
#         job = Job(sim_input=sim_input)
#         # job.go_to_workdir()
#         job.start()
#         try:
#             job_ids, results = func(*args[:], **kwargs)
#             job.complete()
#         except:
#             job.fail()
#             raise
#         job.update_output(results)
#         job.leave_workdir()
#         return job_ids, results
#     return wrapper

# def uses_calc(func):
#     ''' 
#     Step wrapper. This handles things like:
#     - Going to the correct workdir
#     - initializing output files for progress tracking
#     while being minimally invasive
#     '''
#     # functools makes it so that the docs work instead of giving the
#     # wrapper docstring
#     @functools.wraps(func)
#     def wrapper(*args, **kwargs):
#         calc_params = kwargs.get('calc_params', False)
#         assert calc_params, 'Provide calc_params for simulation'
#         print(kwargs)
#         workdir = kwargs.get('workdir', '.')
#         sim_input = {
#             'workdir': workdir,
#         }
#         job = Job(sim_input=sim_input)
#         job.go_to_workdir()
#         job.start()
#         try:
#             job_ids, results = func(*args[:], **kwargs)
#             job.complete()
#         except:
#             job.fail()
#             raise
#         job.update_output(results)
#         job.leave_workdir()
#         return job_ids, results
#     return wrapper

# def uses_calc(func):
#     ''' 
#     Step wrapper. This handles things like:
#     - Going to the correct workdir
#     - initializing output files for progress tracking
#     while being minimally invasive
#     '''
#     # functools makes it so that the docs work instead of giving the
#     # wrapper docstring
#     @functools.wraps(func)
#     def wrapper(*args, **kwargs):
#         env_input = kwargs.get('env_input', {})
#         calc_input = kwargs.get('calc_input', False)
#         sim_input = kwargs.get('sim_input', {})
#         assert calc_input, 'Provide calc_input for simulation'
#         workdir = kwargs.get('workdir', '.')
#         calc_id = kwargs.get('calc_id', False)
#         assert calc_id, 'Provide calc_id in sim_input'
#         calc_params = calc_input[calc_id]
#         job = UnitJob(
#             sim_input=sim_input,
#             env_input=env_input,
#             calc_input=calc_input
#         )
#         job.go_to_workdir()
#         job.start()
#         try:
#             job_ids, results = func(calc_params, *args[:], **kwargs)
#             job.complete()
#         except:
#             job.fail()
#             raise
#         job.update_output(results)
#         job.leave_workdir()
#         return job_ids, results
#     return wrapper

# def branch(func):
#     ''' 
#     Step wrapper. This handles things like:
#     - Going to the correct workdir
#     - initializing output files for progress tracking
#     while being minimally invasive
#     '''
#     # functools makes it so that the docs work instead of giving the
#     # wrapper docstring
#     @functools.wraps(func)
#     def wrapper(*arg, **kwargs):
#         workdir = kwargs.get('workdir', '.')
#         prefix = kwargs.get('prefix', '')
#         sim_input = {
#             'workdir': workdir,
#             'prefix': prefix,
#         }
#         job = Job(*arg, sim_input=sim_input)
#         job.go_to_workdir()
#         job.start()
#         try:
#             job_ids, results = func(*arg, **kwargs)
#         except:
#             job.fail()
#             raise
#         job.update_output(results)
#         job.complete()
#         job.leave_workdir()
#         return job_ids, results
#     return wrapper

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
