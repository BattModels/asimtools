#!/usr/bin/env python
'''
Calculates EOS

Author: mkphuthi@github.com
'''

import importlib
import sys
import getopt
from asimtools.utils import read_yaml

def parse_command_line(argv) -> tuple[dict, dict]:
    ''' Parse command line inputs for simulation script '''
    usage = 'python singlepoint.py calc_params.yaml sim_params.yaml [-h]'

    assert len(argv) >= 2, usage
    try:
        calc_inputfile = argv[0]
        calc_params = read_yaml(calc_inputfile)
    except IndexError:
        print('Failed to load calc_inputfile\n', usage)
        raise

    try:
        sim_inputfile = argv[1]
        sim_params = read_yaml(sim_inputfile)
    except FileNotFoundError:
        print(f'Failed to load {sim_inputfile}\n', usage)
        raise

    try:
        opts, _ = getopt.getopt(
            argv[2:],
            "h",
            []
        )
    except getopt.GetoptError:
        print(usage)
        raise

    for opt, _ in opts:
        if opt == '-h':
            print(usage)
            sys.exit()

    return calc_params, sim_params


def main(argv):
    ''' Main '''
    calc_input, sim_input = parse_command_line(argv)
    script = sim_input['script']
    sim_module = importlib.import_module(script)
    sim_func = getattr(sim_module, script)
    sim_func(
        calc_input,
        **sim_input,
    )


if __name__ == "__main__":
    main(sys.argv[1:])
