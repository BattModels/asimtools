#!/usr/bin/env python
'''
Calculates EOS

Author: mkphuthi@github.com
'''

import importlib
from typing import TypeVar, Tuple
import sys
import argparse
import getopt
from asimtools.utils import read_yaml
from asimtools import scripts

def parse_command_line(args) -> Tuple[dict, dict]:
    parser = argparse.ArgumentParser()
    parser.add_argument('calc_inputfile', metavar='calculator_configuration_file', type=str,
                        help='active_learning configuration file')
    parser.add_argument('sim_inputfile', metavar='simulation_configuration_file', type=str,
                        help='active_learning configuration file')
    args = parser.parse_args()

    calc_params = read_yaml(args.calc_inputfile)
    sim_params = read_yaml(args.sim_inputfile)

    return calc_params, sim_params


def main(args=None):
    ''' Main '''
    calc_input, sim_input = parse_command_line(args)
    script = sim_input['script']
    
    sim_module = importlib.import_module('.'+script,'asimtools.scripts')
    sim_func = getattr(sim_module, script)
    sim_func(
        calc_input,
        **sim_input,
    )


if __name__ == "__main__":
    main(sys.argv[1:])
