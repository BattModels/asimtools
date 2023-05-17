"""
Tests for running scripts using asim_run.py
"""
import os
import yaml

from asimtools.scripts.asim_run import parse_command_line

def test_asim_run_singlepoint():

    calculator_file_name = 'calc.yaml'
    simulation_file_name = 'sim.yaml'

    calc_config = {
        'name': 'LennardJones',
        'job': {
            'slurm': 'true',
            'batch': 'false',
            'interactive': 'true',
        },
    }
    sim_config = {
        'script': 'singlepoint',
        'image': {
            'symbol': 'Ar',
            'crystal_structure': 'fcc',
            'a': 3.4,
        },
        'prefix': 'Ar_LJ',
        'properties': ['energy', 'forces', 'stress'],
    }

    with open(calculator_file_name,'w',encoding='utf8') as fl:
        yaml.dump(calc_config,fl)
    with open(simulation_file_name,'w',encoding='utf8') as fl:
        yaml.dump(sim_config,fl)

    calc_config_parse, sim_config_parse = parse_command_line(
        [calculator_file_name,simulation_file_name]
    )

    ## some checks
    assert calc_config == calc_config_parse
    assert sim_config == sim_config_parse

    ### cleanup
    os.remove(calculator_file_name)
    os.remove(simulation_file_name)
if __name__ == '__main__':
    test_asim_run_singlepoint()
