'''
Set of tools for handling jobs in slurm, interactively or in the teminal for
a specific calculator.

Author: mkphuthi@github.com
'''

import subprocess
import os
import pkgutil
from pathlib import Path
from datetime import datetime
import logging
from glob import glob
from typing import List, TypeVar, Dict, Tuple, Union
from copy import deepcopy
from colorama import Fore
import ase.io
from ase.parallel import paropen
import numpy as np
from asimtools.utils import (
    read_yaml,
    write_yaml,
    join_names,
    get_atoms,
    get_images,
    get_env_input,
    get_calc_input,
    get_logger,
    check_if_slurm_job_is_running,
)

Atoms = TypeVar('Atoms')

START_MESSAGE = '+' * 15 + ' ASIMTOOLS START' + '+' * 15 + '\n'
STOP_MESSAGE = '+' * 15 + ' ASIMTOOLS STOP' + '+' * 15 + '\n'

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
        self.workdir = Path(sim_input.get('workdir', '.')).resolve()

        # sim_input will be in workdir from now on so set workdir as a relative
        # path, but in this class/subclasses, always use self.workdir
        # which is the full path to avoid confusion
        self.sim_input['workdir'] = './'
        self.launchdir = Path(os.getcwd())
        if self.sim_input.get('src_dir', False):
            self.sim_input['src_dir'] = self.launchdir

        self.env_id = self.sim_input.get('env_id', None)
        if self.env_id is not None:
            self.env = self.env_input[self.env_id]
        else:
            self.env = {
                'mode': {
                    'use_slurm': False,
                    'interactive': True,
                }
            } # Default is to run in current console

    def mkworkdir(self) -> None:
        ''' Creates the work directory if it doesn't exist '''
        if not self.workdir.exists():
            self.workdir.mkdir()

    def update_status(self, status: str) -> None:
        ''' Updates job status to specificied value '''
        self.update_output({'status': status})

    def start(self) -> None:
        ''' Updates the output to signal that the job was started '''
        self.update_output({
            'start_time': datetime.now().strftime('%H:%M:%S, %m/%d/%y')
        })
        self.update_status('started')

    def complete(self) -> None:
        ''' Updates the output to signal that the job was started '''
        self.update_output({
            'end_time': datetime.now().strftime('%H:%M:%S, %m/%d/%y')
        })
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

    def get_output_yaml(self) -> Dict:
        ''' Get current output file '''
        out_fname = 'output.yaml'
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

    def update_sim_input(self, new_params) -> None:
        ''' Update simulation parameters '''
        self.sim_input.update(new_params)

    def update_calc_input(self, new_params) -> None:
        ''' Update calculator parameters '''
        self.calc_input.update(new_params)

    def update_env_input(self, new_params) -> None:
        ''' Update calculator parameters '''
        self.env_input.update(new_params)
        self.env_id = self.sim_input.get('env_id', None)
        if self.env_id is not None:
            self.env = self.env_input[self.env_id]

    def set_workdir(self, workdir) -> None:
        ''' Set working directory both in sim_input and instance '''
        workdir = Path(workdir)
        self.sim_input['workdir'] = str(workdir.resolve())
        self.workdir = workdir

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
            print(f'Status: {status}')
        return complete, status

    def update_output(self, output_update: Dict) -> None:
        ''' Update output.yaml if it exists or write a new one '''
        output = self.get_output()
        output.update(output_update)
        write_yaml(self.get_output_yaml(), output)

    def get_logger(self, logfilename='job.log', level='info'):
        ''' Get the logger '''
        assert self.workdir.exists(), 'Work directory does not exist yet'
        logger = get_logger(logfile=self.workdir / logfilename, level=level)
        return logger

    def _gen_slurm_batch_preamble(
        self,
        slurm_params=None,
        extra_flags=None
    ) -> None:
        ''' Generate the txt with job configuration but no run commands '''
        if slurm_params is None:
            slurm_params = self.env.get('slurm', {})

        txt = '#!/usr/bin/sh\n\n'
        flags = slurm_params.get('flags', [])
        if isinstance(flags, dict):
            flag_list = []
            for k, v in flags.items():
                assert k.startswith('-'), \
                    'Slurm flags must start with "-" or "--"'
                if k.startswith('--'):
                    delim = '='
                else:
                    delim = ' '

                flag_list.append(delim.join([k, str(v)]))

            flags = flag_list

        jobname_flag = '-J ' + str(self.sim_input.get('asimmodule', 'job'))
        for i, flag in enumerate(flags):
            if flag.strip().startswith('-J') or flag.strip().startswith('--job-name'):
                jobname_flag = flags.pop(i)
                break

        # Allow naming jobs in sim_input files
        job_name = self.sim_input.get('job_name', False)
        if job_name:
            jobname_flag = f'-J {job_name}'

        flags.append(jobname_flag)
        if extra_flags is not None:
            # Only supporting list for extraflags, it is an internal variable
            flags.extend(extra_flags)
        for flag in flags:
            txt += f'#SBATCH {flag}\n'

        return txt

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

        # Check if the asimmodule being called uses a calc_id to
        # get precommands, postcommands, run_prefixes and run_suffixes
        self.calc_id = self.sim_input.get('args', {}).get('calc_id', None)
        if self.calc_id is not None:
            self.calc_params = self.calc_input[self.calc_id]
        else:
            self.calc_params = {}

    def gen_run_command(self) -> str:
        '''
        Generates the command to run the job, independent of slurm
        '''
        asimmodule = self.sim_input.get('asimmodule', False)
        assert asimmodule, 'Specify a asimmodule in sim_input'

        mode_params = self.env.get('mode', {
            'use_slurm': False, 'interactive': True
        })
        env_run_prefix = mode_params.get('run_prefix', '')
        env_run_suffix = mode_params.get('run_suffix', '')

        calc_run_prefix = self.calc_params.get('run_prefix', '')
        calc_run_suffix = self.calc_params.get('run_suffix', '')

        if self.sim_input.get('debug', False):
            debug_flag = ' -d'
        else:
            debug_flag = ''

        # GPAW requires python asimmodule to be called with gpaw python
        # For now this is a work around but perhaps there exists a cleaner way
        if self.calc_params.get('name', '') == 'GPAW':
            asloader = pkgutil.get_loader('asimtools.scripts.asim_run')
            asfile = asloader.get_filename()
            alias = f'gpaw python {asfile}'
        else:
            alias = 'asim-run'

        txt = f'{env_run_prefix} {calc_run_prefix} '
        txt += f'{alias} sim_input.yaml --calc=calc_input.yaml' + debug_flag
        txt += f' {env_run_suffix} {calc_run_suffix} '
        return txt

    def _gen_slurm_script(self, write: bool = True) -> None:
        '''
        Generates a slurm job file and if necessary writes it in 
        the work directory for the job
        '''
        slurm_params = self.env.get('slurm', {})
        calc_params = self.calc_params

        txt = self._gen_slurm_batch_preamble()
        txt += '\n'
        for command in slurm_params.get('precommands', []):
            txt += command + '\n'
        for command in calc_params.get('precommands', []):
            txt += command + '\n'
        txt += '\n'
        txt += 'echo "Job started on `hostname` at `date`"\n'
        txt += self.gen_run_command() + '\n'
        for command in slurm_params.get('postcommands', []):
            txt += command + '\n'
        for command in calc_params.get('postcommands', []):
            txt += command + '\n'
        txt += 'echo "Job ended at `date`"'

        if write:
            slurm_file = self.workdir / 'job.sh'
            slurm_file.write_text(txt)

    def _gen_slurm_interactive_txt(self) -> str:
        '''
        Generates an srun command to run the job
        '''
        slurm_params = self.env.get('slurm', False)

        txt = 'srun --unbuffered ' # unbuffered to avoid IO lag
        txt += ' '.join(slurm_params.get('flags'))
        txt += self.gen_run_command()

        return txt

    def gen_input_files(
        self,
        write_calc_input: bool = True,
        write_env_input: bool = True,
        # use_links: bool = True,
        write_image: bool = True,
    ) -> None:
        ''' Write input files to working directory '''
        workdir = self.workdir

        # This is the sim_input that will be written to disk so it will have
        # slightly different paths and the image(s) will be image_files
        sim_input = deepcopy(self.sim_input)
        # The sim_input file will already be in the work directory
        sim_input['workdir'] = '.'
        # Collect useful variables
        mode_params = self.env.get('mode', {})
        use_slurm = mode_params.get('use_slurm', False)
        interactive = mode_params.get('interactive', True)
        overwrite = self.sim_input.get('overwrite', False)
        output_file = workdir / 'output.yaml'
        sim_input_file = workdir / 'sim_input.yaml'

        # Check if the job already exists or not and if to overwrite
        if not workdir.exists():
            workdir.mkdir()
        elif output_file.exists() and sim_input_file.exists():
            complete = read_yaml(output_file).get('end_time', False)
            if complete and not overwrite:
                txt = f'Skipped writing in {self.workdir}, files/directory '
                txt += 'already exist. Set overwrite=True to overwrite files '
                txt += 'or delete results directory first (recommended)'
                self.get_logger().warning(txt)
                print(Fore.RED + 'WARNING: ' + txt + Fore.WHITE)
                return None

        # If the sim_input uses an image/images, we want to write it in
        # the workdir and point to it using a filename. This is so that
        # there is a record of the input image for debugging and the input
        # method is consistent with complex workflows
        image = self.sim_input.get('args', {}).get('image', False)

        # write_image is to account for chains where the image may not have
        # been produced yet by a previous step, otherwise always read and write
        if image and write_image:
            atoms = get_atoms(**image)
            input_image_file = 'image_input.xyz' # Relative to workdir
            # input_image_file = self.workdir.resolve() / 'image_input.xyz'
            atoms.write(
                self.workdir / input_image_file,
                format='extxyz'
            )
            sim_input['args']['image'] = {
                'image_file': str(input_image_file),
                # 'image_file': str((self.workdir / input_image_file).resolve())
            }

        images = self.sim_input.get('args', {}).get('images', False)
        if images and write_image:
            if images.get('images', False) or images.get('image_file', False):
                images = get_images(**images)
                input_images_file = 'images_input.xyz' # Relative to workdir
                ase.io.write(
                    self.workdir / input_images_file,
                    images,
                    format='extxyz',
                )
                sim_input['args']['images'] = {
                    'image_file': str(input_images_file),
                }

        write_yaml(workdir / 'sim_input.yaml', sim_input)
        if write_calc_input:
            write_yaml(workdir / 'calc_input.yaml', self.calc_input)
        if write_env_input:
            write_yaml(workdir / 'env_input.yaml', self.env_input)
        write_yaml(self.get_output_yaml(), {'status': 'clean'})

        if use_slurm and not interactive:
            self._gen_slurm_script()
        return None

    def submit(
        self,
        dependency: Union[List,None] = None,
        write_image: bool = True,
    ) -> Union[None,List[str]]:
        '''
        Submit a job using slurm, interactively or in the terminal
        '''

        msg = START_MESSAGE
        msg += f'asimmodule: {self.sim_input["asimmodule"]}\n'
        msg += f'Work directory: {self.workdir}'
        logmsg = f'Submitting "{self.sim_input["asimmodule"]}" in "{self.workdir}"'
        logging.info(logmsg)
        print(msg) # For printing to console

        self.gen_input_files(write_image=write_image)

        logger = self.get_logger()

        _, status = self.get_status(display=True)
        if status in ('discard', 'complete'):
            txt = f'Current job status in "{self.workdir}" is "{status}", '
            txt += 'skipping submission, set "overwrite=True" to ignore'
            logger.warning(txt)
            if status in ('complete'):
                print('Job already completed, not submitting again')
                print(STOP_MESSAGE)
            return None

        if not self.sim_input.get('submit', True):
            logger.warning('Keyword submit=False, skipping submission in \
                %s', self.workdir)
            return None

        cur_dir = Path('.').resolve()
        os.chdir(self.workdir)
        mode_params = self.env.get('mode', {})
        sim_input = self.sim_input
        use_slurm = mode_params.get('use_slurm', False)
        interactive = mode_params.get('interactive', True)

        # Case when running batch jobs, handles job dependency
        if use_slurm and not interactive:
            if dependency is not None:
                dependstr = None
                for dep in dependency:
                    if dependstr is None:
                        dependstr = str(dep) if dep is not None else None
                    else:
                        dependstr += f':{dep}' if dep is not None else ''

                command = ['sbatch', '-d', f'afterok:{dependstr}']
                command += ['--kill-on-invalid-dep=yes', 'job.sh']
            else:
                command = ['sbatch', 'job.sh']

        # Interactive jobs launched with srun
        elif use_slurm and interactive:
            command = self._gen_slurm_interactive_txt().split()

        # Launching the job inline
        else:
            command = self.gen_run_command().split()

        run_job = False
        overwrite = sim_input.get('overwrite', False)
        _, status = self.get_status()
        if overwrite or status == 'clean':
            run_job = True

        if run_job:
            completed_process = subprocess.run(
                command, check=False, capture_output=True, text=True,
            )

            with paropen('stdout.txt', 'w', encoding='utf-8') as output_file:
                output_file.write(completed_process.stdout)

            if completed_process.stderr is not None:
                with paropen('stderr.txt', 'w', encoding='utf-8') as err_file:
                    err_file.write(completed_process.stderr)

            if completed_process.returncode != 0:
                err_msg = f'See {self.workdir / "stderr.txt"} for traceback.'
                logger.error(err_msg)
                completed_process.check_returncode()

        os.chdir(cur_dir)

        # Return job ids for chaining dependencies in batch mode
        if use_slurm and not interactive and run_job:
            job_ids = [int(completed_process.stdout.split(' ')[-1])]
            self.update_output({'job_ids': job_ids})
            return job_ids

        logmsg = f'Submitted "{self.sim_input["asimmodule"]}" in "{self.workdir}"'
        msg = 'Successfully submitted\n'
        msg += STOP_MESSAGE
        logger.info(logmsg)
        print(msg)

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
        # with "id-{sim_id}" followed by user labels
        new_sim_input = {}
        njobs = len(sim_input)
        num_digits = 4
        assert njobs < 1000, \
            f'ASIMTools id_nums are limited to {num_digits} digits, found \
            {njobs} jobs! Having that many jobs is not very storage efficient.'

        sim_id_changed = False
        for i, (sim_id, subsim_input) in enumerate(sim_input.items()):
            assert 'asimmodule' in subsim_input, 'Check sim_input format,\
                must have asimmodule and each job with a unique key'

            # We set a fixed number of digits so that the slurm array ID
            # matches the directory name for easy resubmission of arrays
            prefix = 'id-'+f'{i}'.rjust(num_digits, '0')
            if not sim_id.startswith(prefix):
                sim_id_changed = True

        if sim_id_changed:
            for i, (sim_id, subsim_input) in enumerate(sim_input.items()):
                prefix = 'id-'+f'{i}'.rjust(num_digits, '0')
                new_sim_id = join_names([prefix,sim_id])
                new_subsim_input = deepcopy(sim_input[sim_id])
                new_subsim_input['distributed_id'] = i
                new_sim_input[new_sim_id] = new_subsim_input
            sim_input = new_sim_input
        unitjobs = []
        for sim_id, subsim_input in sim_input.items():
            assert 'asimmodule' in subsim_input, 'Check sim_input format,\
                must have asimmodule and each job with a unique key'

            subsim_input['workdir'] = sim_id
            unitjob = UnitJob(subsim_input, env_input, calc_input)
            unitjobs.append(unitjob)

        # If all the jobs have the same config and use slurm, use a job array
        env_id = unitjobs[0].env_id
        same_env = np.all(
            [(uj.env_id == env_id) for uj in unitjobs]
        )

        all_slurm = np.all(
            [uj.env['mode'].get('use_slurm', False) for uj in unitjobs]
        )

        if same_env and all_slurm:
            self.use_array = True
        else:
            self.use_array = False

        self.unitjobs = unitjobs

    def submit_jobs(
        self,
        **kwargs, # Necessary for compatibility with job.submit
    ) -> Union[None,List[str]]:
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
        **kwargs, # Necessary for compatibility with job.submit
    ) -> Union[None,List[str]]:
        '''
        Submits a job array if all the jobs have the same env and use slurm
        '''
        self.gen_input_files()
        logger = self.get_logger()

        unitjobs = [] # Track jobs that are supposed to be submitted
        for unitjob in self.unitjobs:
            if unitjob.sim_input.get('submit', True):
                unitjobs.append(unitjob)

        njobs = len(unitjobs)
        if njobs == 0:
            return None

        # For limiting number of jobs launched in array
        if array_max is not None:
            arr_max_str = f'%{array_max}'
        else:
            arr_max = unitjobs[0].env['mode'].get('array_max', False)
            if arr_max:
                arr_max_str = f'%{arr_max}'
            else:
                arr_max_str = ''

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

        with paropen('stdout.txt', 'w', encoding='utf-8') as output_file:
            output_file.write(completed_process.stdout)

        if completed_process.stderr is not None:
            with paropen('stderr.txt', 'w', encoding='utf-8') as err_file:
                    err_file.write(completed_process.stderr)

        if completed_process.returncode != 0:
            err_msg = f'See {self.workdir / "stderr.txt"} for traceback.'
            logger.error(err_msg)
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

        txt = self._gen_slurm_batch_preamble(
            slurm_params=slurm_params,
            extra_flags=[
                '-o slurm_stdout.id-%a_j%A',
                '-e slurm_stderr.id-%a_j%A',
            ]
        )

        txt += '\nif [[ ! -z ${SLURM_ARRAY_TASK_ID} ]]; then\n'
        txt += '    fls=( id-* )\n'
        txt += '    WORKDIR=${fls[${SLURM_ARRAY_TASK_ID}]}\n'
        txt += 'fi\n\n'
        txt += 'CUR_DIR=`pwd`\n'
        txt += 'cd ${WORKDIR}\n'
        txt += '\n'
        txt += '\n'.join(slurm_params.get('precommands', []))
        txt += '\n'
        txt += '\n'.join(self.unitjobs[0].calc_params.get('precommands', []))
        txt += '\n'
        txt += 'echo "LAUNCHDIR: ${CUR_DIR}"\n'
        txt += 'echo "WORKDIR: ${WORKDIR}"\n'
        txt += 'echo "Job started on `hostname` at `date`"\n'
        txt += self.unitjobs[0].gen_run_command() + '\n'
        txt += '\n'.join(slurm_params.get('postcommands', []))
        txt += '\n'
        txt += 'cd ${CUR_DIR}\n'
        txt += 'echo "Job ended at `date`"'

        if write:
            slurm_file = self.workdir / 'job_array.sh'
            slurm_file.write_text(txt)

    def gen_input_files(self) -> None:
        ''' Write input files to working directory '''

        if self.use_array:
            self._gen_array_script()
        for unitjob in self.unitjobs:
            unitjob.gen_input_files()

    def submit(self, **kwargs) -> None:
        ''' 
        Submit a job using slurm, interactively or in the terminal
        '''

        cur_dir = Path('.').resolve()
        os.chdir(self.workdir)
        job_ids = None
        if self.use_array:
            job_ids = self.submit_array(**kwargs)
        else:
            job_ids = self.submit_jobs(**kwargs)
        self.update_output({'job_ids': job_ids}) # For chained job
        os.chdir(cur_dir)
        return job_ids

class ChainedJob(Job):
    '''
    Jobs made of smaller sequential unitjobs. Note that this only works
    well with unit jobs. If one of the UnitJobs calls asimmodules, internally,
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

        for key in sim_input:
            assert key.startswith('step-'), \
                f'Keys take the form "step-[step_id]" from id=0 not "{key}"'

        # Each key in sim_input should be in the form "step-id" where id
        # is replaced by the order in which the commands are run
        unitjobs = []
        for i in range(len(sim_input)):
            step_id = f'step-{i}'
            subsim_input = deepcopy(sim_input[step_id])
            unitjob = UnitJob(subsim_input, env_input, calc_input)
            unitjob.set_workdir(step_id)
            unitjobs.append(unitjob)

        self.unitjobs = unitjobs

    def get_last_output(self) -> Dict:
        ''' Returns the output of the last job in the chain '''
        return self.unitjobs[-1].get_output()

    def submit(self, dependency: Union[List,None] = None) -> List:
        ''' 
        Submit a job using slurm, interactively or in the terminal
        '''
        cur_dir = Path('.').resolve()
        os.chdir(self.workdir)
        logger = self.get_logger()
        step = 0 #self.get_current_step() TODO: This feature is not used yet
        are_interactive_jobs = [
            uj.env['mode'].get('interactive', False) \
                for uj in self.unitjobs
        ]
        all_interactive = np.all(are_interactive_jobs)
        all_not_interactive = np.sum(are_interactive_jobs) == 0
        assert all_interactive or all_not_interactive, \
            'Jobs in chain must  all be interactive=True or interactive=False'

        are_using_slurm = [
            uj.env['mode'].get('use_slurm', False) \
                for uj in self.unitjobs
        ]
        all_slurm = np.all(are_using_slurm)
        all_not_slurm = np.sum(are_using_slurm) == 0

        assert all_slurm or all_not_slurm, \
            'Jobs in chain must all have use_slurm=True or use_slurm=False'

        # Only use slurm dependencies if all the jobs are batch jobs
        if all_slurm and all_not_interactive:
            if step != len(self.unitjobs):
                dependency = None
                for i, unitjob in enumerate(self.unitjobs[step:]):
                    cur_step_complete, status = unitjob.get_status()
                    if not cur_step_complete:
                        slurm_flags = unitjob.env['slurm']['flags']
                        if isinstance(slurm_flags, list):
                            job_name = unitjob.sim_input.get('job_name',False)
                            if not job_name:
                                job_name = f'-J step-{step+i}'
                            else:
                                job_name = f'-J step-{step+i}-{job_name}'
                            unitjob.env['slurm']['flags'].append(
                                job_name
                            )
                        elif isinstance(slurm_flags, dict):
                            unitjob.env['slurm']['flags']['-J'] = \
                                f'step-{step+i}'

                        write_image = False
                        # Write image first step in chain being run/continued
                        if i == 0:
                            write_image = True
                        dependency = unitjob.submit(
                            dependency=dependency,
                            write_image=write_image,
                        )
                    else:
                        txt = f'Skipping step-{step+i} since '
                        txt += f'status is {status}'
                        logger.warning(txt)

            job_ids = dependency

        # Otherwise just submit them one after the other
        # We only write the image if it's the first job, otherwise we refer to 
        # wherever the image comes from in case it has to come from a previous
        # step
        else:
            for i, unitjob in enumerate(self.unitjobs[step:]):
                if i == 0:
                    write_image = True
                else:
                    write_image = False

                job_ids = unitjob.submit(write_image=write_image)

        os.chdir(cur_dir)

        return job_ids


def load_job_from_directory(workdir: str):
    ''' Loads a job from a given directory '''
    workdir = Path(workdir)
    assert workdir.exists(), f'Work director {workdir} does not exist'
    logger = get_logger()
    sim_inputs = glob(str(workdir / 'sim_input.yaml'))
    if len(sim_inputs) != 1:
        logger.error('Multiple or no sim_input.yaml files in %s', {workdir})
    try:
        sim_input = read_yaml(glob(str(workdir / 'sim_input.yaml'))[0])
    except IndexError as exc:
        logger.error('sim_input.yaml not found in %s', {str(workdir)})
        raise exc

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

    # This makes sure that wherever we may be loading the job from, we refer
    # to the correct input/output files. As of now, it does not seem necessary
    # to also change the workdir in job.sim_input for consistency
    job.workdir = Path(workdir)
    return job

def create_unitjob(sim_input, env_input, workdir, calc_input=None):
    """Helper for making a generic UnitJob object, mostly for testing"""
    env_id = list(env_input.keys())[0]
    sim_input['env_id'] = env_id
    if calc_input is not None:
        calc_id = list(calc_input.keys())[0]
        sim_input['calc_id'] = calc_id
    sim_input['workdir'] = workdir
    unitjob = UnitJob(
        sim_input,
        env_input,
        calc_input
    )
    return unitjob
