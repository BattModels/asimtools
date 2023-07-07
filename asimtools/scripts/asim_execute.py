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
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument(
        'sim_input_file',
        metavar='simulation_input_file',
        type=str,
        help='simulation input yaml file',
    )
    parser.add_argument(
        '-e',
        '--env',
        metavar='environment_input_file',
        type=str,
        help='environment input yaml file',
    )
    parser.add_argument(
        '-c',
        '--calc',
        metavar='calculator_input_file',
        type=str,
        help='calculator input yaml file',
    )
    args = parser.parse_args(args)

    sim_input = read_yaml(args.sim_input_file)
    calc_input = args.calc
    env_input = args.env
    if env_input is not None:
        env_input = read_yaml(env_input)

    calc_input = args.calc
    if calc_input is not None:
        calc_input = read_yaml(calc_input)

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
    # job.start()
    try:
        job.submit()
    except:
        job.fail()
        raise
    # job.update_output({'job_ids': job_ids})
    job.leave_workdir()


if __name__ == "__main__":
    main(sys.argv[1:])
