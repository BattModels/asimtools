#!/usr/bin/env python
'''
Execute a asimmodule given the sim_input.yaml and optionally,
a calc_input.yaml. The called asimmodule will be run directly in the current
directory and environment
'''

import importlib
import sys
import os
from pathlib import Path
import argparse
import subprocess
from typing import Dict, Tuple
from asimtools.utils import read_yaml, get_logger
from asimtools.job import load_job_from_directory


def parse_command_line(args) -> Tuple[Dict, str]:
    ''' Parse command line input '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'sim_input_file',
        metavar='simulation_input_file',
        type=str,
        help='simulation input yaml file'
    )
    parser.add_argument(
        '-c',
        '--calc',
        metavar='calculator_input_file',
        type=str,
        help='calculator input yaml file'
    )
    parser.add_argument(
        '-d',
        '--debug',
        action='store_true',
        help='Set logging level'
    )
    args = parser.parse_args(args)

    sim_input = read_yaml(args.sim_input_file)
    if args.debug:
        sim_input['debug'] = True
    calc_input_file = args.calc

    return sim_input, calc_input_file


def main(args=None) -> None:
    ''' Main '''
    sim_input, calc_input_file = parse_command_line(args)

    if sim_input.get('debug', False):
        level = 'debug'
    else:
        level = 'warning'
    logger = get_logger(level=level)
    logger.debug('Entered asim_run')

    old_calc_input_var = os.getenv("ASIMTOOLS_CALC_INPUT", None)
    if calc_input_file is not None:
        assert Path(calc_input_file).exists(), 'Specify valid calc input file'
        os.environ["ASIMTOOLS_CALC_INPUT"] = calc_input_file
    asimmodule = sim_input['asimmodule']
    precommands = sim_input.get('precommands', [])

    for precommand in precommands:
        command = precommand.split()
        completed_process = subprocess.run(command, check=True)
        completed_process = subprocess.run(
            command, check=False, capture_output=True, text=True,
        )

        if completed_process.returncode != 0:
            logger.error('asim-run: Failed to run precommand %s', precommand)
            logger.error(completed_process.stderr)
            err_txt = f'ERROR: sim_input precommand "{" ".join(command)}" '
            err_txt += 'failed as above. Make sure precommand can be '
            err_txt += "executed by python's subprocess.run() routine"
            logger.error(err_txt)

            completed_process.check_returncode()

    module_name = asimmodule.split('/')[-1].split('.py')[0]
    func_name = module_name.split('.')[-1] # must match asimmodule name

    # Check if the asimmodule is a core asimmodule
    spec = importlib.util.find_spec('.'+module_name,'asimtools.asimmodules')
    if spec is not None:
        sim_module = importlib.import_module(
            '.'+module_name,'asimtools.asimmodules'
        )
        logger.debug('Loaded core asimmodule %s', module_name)
    else:
        # Try and find a asimmodule from the set asimmodule directory
        asimmodule_dir = os.environ.get('ASIMTOOLS_ASIMMODULE_DIR', False)
        if asimmodule_dir:
            # asimmodule = asimmodule.replace('.', '/')
            asimmodule = Path(asimmodule_dir) / asimmodule
            spec = importlib.util.spec_from_file_location(
                    module_name, asimmodule
                )
            logger.debug('Loading asimmodule from asimmodule directory %s', asimmodule)
        if spec is None:
            spec = importlib.util.spec_from_file_location(
                module_name, asimmodule
            )
            logger.debug('Loading asimmodule from full path %s', asimmodule)

        try:
            sim_module = importlib.util.module_from_spec(spec)
        except Exception as exc:
            err_msg = f'Failed to load module for asimmodule: "{asimmodule}"'
            logger.error(err_msg)
            raise AttributeError(err_msg) from exc

        sys.modules[module_name] = sim_module
        try:
            spec.loader.exec_module(sim_module)
        except Exception as exc:
            t = f'Failed to load asimmodule "{asimmodule}". Possible causes:\n'
            t += '* The asimmodule has a bug, see traceback above\n'
            t += '* ASIMTOOLS_ASIMMODULE_DIR variable not set properly.\n'
            t += '* You can provide the full path to the asimmodule'
            logger.error(t)
            raise FileNotFoundError(t) from exc

    sim_func = getattr(sim_module, func_name)

    cwd = Path('.').resolve()
    job = load_job_from_directory(cwd)
    job.start()

    try:
        results = sim_func(**sim_input.get('args', {}))
        logger.info('Successfully ran asimmodule "%s" in "%s"', asimmodule, cwd)
    except:
        logger.info('Failed to run asimmodule "%s" in "%s"', asimmodule, cwd)
        # Restore the CALC_INPUT environment variable
        if old_calc_input_var is not None:
            os.environ["ASIMTOOLS_CALC_INPUT"] = old_calc_input_var
        else:
            if os.getenv("ASIMTOOLS_CALC_INPUT", None) is not None:
                del os.environ["ASIMTOOLS_CALC_INPUT"]
        # if world.rank == 0:
        job.fail()
        raise

    # Restore the previous calc_input in case other srcipts depend on it
    # This is pretty messy right now, could probably be better
    if old_calc_input_var is not None:
        os.environ["ASIMTOOLS_CALC_INPUT"] = old_calc_input_var
    else:
        if os.getenv("ASIMTOOLS_CALC_INPUT", None) is not None:
            del os.environ["ASIMTOOLS_CALC_INPUT"]

    # Get job IDs from internally run asimmodules if any
    job_ids = results.get('job_ids', False)
    # Otherwise get job_ids from wherever this asimmodule is running
    if not job_ids:
        job_ids = os.getenv('SLURM_JOB_ID')
        results['job_ids'] = job_ids

    job.update_output(results)
    job.complete()

    postcommands = sim_input.get('postcommands', [])
    for postcommand in postcommands:
        command = postcommand.split()
        completed_process = subprocess.run(command, check=True)

    logger.debug('Complete')

if __name__ == "__main__":
    main(sys.argv[1:])
