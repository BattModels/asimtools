#!/usr/bin/env python
'''
Execute a workflow given the sim_input.yaml and optionally,
a calc_input.yaml and/or env_input.yaml. The called asimmodule will be run
in the specified workdir and env_id

Author: mkphuthi@github.com
'''

import sys
import shutil
import argparse
from pathlib import Path
from typing import Dict, Tuple
from asimtools.utils import read_yaml
from asimtools.job import UnitJob


def parse_command_line(args) -> Tuple[Dict, Dict, Dict]:
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
    parser.add_argument(
        '-d',
        '--debug',
        action='store_true',
        help='Set logging level to debug'
    )
    parser.add_argument(
        '-f',
        '--force',
        action='store_true',
        help='Remove the work directory before rerunning. Be careful!'
    )
    args = parser.parse_args(args)
    sim_input = read_yaml(args.sim_input_file)

    # Do some cleanup of arguments
    if args.debug:
        sim_input['debug'] = True

    calc_input = args.calc
    env_input = args.env
    if env_input is not None:
        env_input = read_yaml(env_input)

    calc_input = args.calc
    if calc_input is not None:
        calc_input = read_yaml(calc_input)

    workdir = Path(sim_input.get('workdir', 'results'))

    if workdir.exists() and args.force:
        if workdir.resolve() != Path('./').resolve():
            shutil.rmtree(workdir)

    return sim_input, env_input, calc_input

def main(args=None) -> None:
    """Execute a workflow given the sim_input.yaml and optionally,
    a calc_input.yaml and/or env_input.yaml. The called asimmodule will be run
    in the specified workdir and env_id """
    sim_input, env_input, calc_input = parse_command_line(args)
    sim_input['workdir'] = sim_input.get('workdir', 'results')

    job = UnitJob(
        sim_input=sim_input,
        env_input=env_input,
        calc_input=calc_input
    )

    job.gen_input_files(write_calc_input=True, write_env_input=True)
    job.go_to_workdir()
    try:
        job.submit()
    except:
        job.fail()
        raise
    job.leave_workdir()


if __name__ == "__main__":
    main(sys.argv[1:])
