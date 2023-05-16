#!/usr/bin/env python
'''
Calculates EOS

Author: mkphuthi@github.com
'''

import importlib
from typing import Tuple
import sys
import argparse
from asimtools.utils import read_yaml

def parse_command_line(args) -> Tuple[dict, dict]:
    ''' Parse command line input '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'calc_input_file',
        metavar='calculator_configuration_file',
        type=str,
        help='calculator configuration yaml file'
    )
    parser.add_argument(
        'sim_input_file',
        metavar='simulation_configuration_file',
        type=str,
        help='calculator configuration yaml file'
    )
    args = parser.parse_args()

    calc_params = read_yaml(args.calc_input_file)
    sim_params = read_yaml(args.sim_input_file)

    return calc_params, sim_params


def main(args=None):
    ''' Main '''
    calc_input, sim_input = parse_command_line(args)
    script = sim_input['script']
    module_name = script.split('/')[-1].split('.py')[0]
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

    sim_func = getattr(sim_module, module_name)
    sim_func(
        calc_input,
        **sim_input,
    )


if __name__ == "__main__":
    main(sys.argv[1:])
