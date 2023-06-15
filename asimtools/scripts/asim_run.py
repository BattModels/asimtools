#!/usr/bin/env python
'''
Run a script

Author: mkphuthi@github.com
'''

import importlib
import sys
import os
import argparse
from typing import Dict, Tuple
from asimtools.utils import read_yaml, get_calc_input
from asimtools.job import load_job_from_directory


def parse_command_line(args) -> Tuple[Dict, Dict]:
    ''' Parse command line input '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'sim_input_file',
        metavar='simulation_configuration_file',
        type=str,
        help='calculator configuration yaml file'
    )
    parser.add_argument(
        'calc_input_file',
        metavar='calculator_configuration_file',
        type=str,
        help='calculator configuration yaml file'
    )
    args = parser.parse_args(args)

    sim_input = read_yaml(args.sim_input_file)
    if args.calc_input_file != 'none':
        calc_input = read_yaml(args.calc_input_file)
    else:
        calc_input = get_calc_input()

    return sim_input, calc_input


def main(args=None) -> None:
    ''' Main '''
    sim_input, _ = parse_command_line(args)
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
    try:
        results = sim_func(**sim_input['args'])
        success = True
    except:
        print(f'ERROR: Failed to run {func_name}')
        raise
    job = load_job_from_directory('.')
    if success:
        job_id = os.getenv('SLURM_JOB_ID')
        results['job_id'] = job_id
        job.update_output(results)
        job.complete()
    else:
        job.fail()


if __name__ == "__main__":
    main(sys.argv[1:])
