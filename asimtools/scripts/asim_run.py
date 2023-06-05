#!/usr/bin/env python
'''
Calculates EOS

Author: mkphuthi@github.com
'''

import importlib
import sys
from asimtools.utils import parse_command_line

def main(args=None) -> None:
    ''' Main '''
    calc_input, sim_input = parse_command_line(args)
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
    sim_func(
        calc_input,
        **sim_input['args'],
        config_id=sim_input['config_id']
    )


if __name__ == "__main__":
    main(sys.argv[1:])
