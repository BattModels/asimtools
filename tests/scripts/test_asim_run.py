"""
Tests for running scripts using asim_run.py
"""
import os
import yaml

from asimtools.scripts.asim_run import main as asim_run

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

    asim_run([calculator_file_name,simulation_file_name])

    ## some checks
    output_file_name = sim_config['prefix'] + '_output.yaml'
    with open(output_file_name,'r',encoding='utf8') as fl:
        output = yaml.load(fl,Loader=yaml.FullLoader)

    assert output.get('status') == 'complete'

    ### cleanup
    output_image_file_name = sim_config['prefix'] + '_image_output.xyz'
    os.remove(calculator_file_name)
    os.remove(simulation_file_name)
    os.remove(output_file_name)
    os.remove(output_image_file_name)

if __name__ == '__main__':
    test_asim_run_singlepoint()
