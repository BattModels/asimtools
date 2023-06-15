#!/usr/bin/env python
'''
Run a script

Author: mkphuthi@github.com
'''

import sys
import argparse
from typing import Dict, Tuple
from asimtools.utils import read_yaml
from asimtools.job import UnitJob


def parse_command_line(args) -> Tuple[Dict, Dict]:
    ''' Parse command line input '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'sim_input_file',
        metavar='simulation_configuration_file',
        type=str,
        help='calculator configuration yaml file',
    )
    parser.add_argument(
        'env_input_file',
        metavar='environment_configuration_file',
        type=str,
        help='environment configuration yaml file',
    )
    parser.add_argument(
        'calc_input_file',
        metavar='calculator_configuration_file',
        type=str,
        help='calculator configuration yaml file',
    )
    args = parser.parse_args(args)
    #TODO: Take inputs as optional arguments rather than using none
    sim_input = read_yaml(args.sim_input_file)
    if args.env_input_file != 'none':
        env_input = read_yaml(args.env_input_file)
    else:
        env_input = None
    if args.calc_input_file != 'none':
        calc_input = read_yaml(args.calc_input_file)
    else:
        calc_input = None

    return sim_input, env_input, calc_input


def main(args=None) -> None:
    ''' Main '''
    sim_input, env_input, calc_input = parse_command_line(args)
    job = UnitJob(
        sim_input=sim_input,
        env_input=env_input,
        calc_input=calc_input
    )
    job.gen_input_files()
    job.go_to_workdir()
    job.start()
    try:
        job_ids = job.submit()
    except:
        job.fail()
        raise
    # job.update_output({'job_ids': job_ids})
    job.leave_workdir()


if __name__ == "__main__":
    main(sys.argv[1:])
