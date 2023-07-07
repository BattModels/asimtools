#!/usr/bin/env python
'''
Execute a workflow given the sim_input.yaml and optionally,
a calc_input.yaml. The called script will be run directly in the current
directory and environment
'''

import importlib
import sys
import os
from pathlib import Path
import argparse
from typing import Dict, Tuple
from asimtools.utils import read_yaml
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
    args = parser.parse_args(args)

    sim_input = read_yaml(args.sim_input_file)
    calc_input_file = args.calc

    return sim_input, calc_input_file


def main(args=None) -> None:
    ''' Main '''
    sim_input, calc_input_file = parse_command_line(args)
    old_calc_input_var = os.getenv("ASIMTOOLS_CALC_INPUT", None)
    if calc_input_file is not None:
        assert Path(calc_input_file).exists(), 'Specify valid calc input file'
        os.environ["ASIMTOOLS_CALC_INPUT"] = calc_input_file
    script = sim_input['script']
    module_name = script.split('/')[-1].split('.py')[0]
    func_name = module_name.split('.')[-1]
    spec = importlib.util.find_spec('.'+module_name,'asimtools.scripts')

    if spec is not None:
        sim_module = importlib.import_module(
            '.'+module_name,'asimtools.scripts'
        )
    else:
        spec = importlib.util.spec_from_file_location(
            module_name, script
        )
        sim_module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = sim_module
        spec.loader.exec_module(sim_module)

    sim_func = getattr(sim_module, func_name)

    print(f'asim-run: running {func_name}')
    job = load_job_from_directory('.')
    job.start()
    try:
        results = sim_func(**sim_input['args'])
    except:
        print(f'ERROR: Failed to run {func_name}')
        if old_calc_input_var is not None:
            os.environ["ASIMTOOLS_CALC_INPUT"] = old_calc_input_var
        else:
            if os.getenv("ASIMTOOLS_CALC_INPUT", False):
                del os.environ["ASIMTOOLS_CALC_INPUT"]
        job.fail()
        raise

    # Restore the previous calc_input in case other srcipts depend on it
    # This is pretty messy right now, could probably be better
    if old_calc_input_var is not None:
        os.environ["ASIMTOOLS_CALC_INPUT"] = old_calc_input_var
    else:
        if os.getenv("ASIMTOOLS_CALC_INPUT", False):
            del os.environ["ASIMTOOLS_CALC_INPUT"]

    # Get job IDs from internally run scripts if any
    job_ids = results.get('job_ids', False)
    # Otherwise get job_ids from wherever this script is running
    if not job_ids:
        job_ids = os.getenv('SLURM_JOB_ID')
        results['job_ids'] = job_ids
    job.update_output(results)
    job.complete()

if __name__ == "__main__":
    main(sys.argv[1:])
