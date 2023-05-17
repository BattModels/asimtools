"""
Tests for running scripts using asim_run.py
"""
import os
from asimtools.utils import write_yaml, read_yaml
from asimtools.scripts.asim_run import main as asim_run

# Path of script that does nothing
SCRIPT_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    '../data/scripts/do_nothing_script.py',
)

def test_asim_run_singlepoint(
    tmp_path,
    lj_interactive_calc_input,
    singlepoint_sim_input
):
    ''' 
    Test asim_run ability to launch and complete a script using LJ
    '''
    # Create input files in temporary directory that will 
    # be automatically deleted using the fixtures in conftest.py
    os.chdir(tmp_path)
    calc_input_file = 'calc_input.yaml'
    sim_input_file = 'sim_input.yaml'
    write_yaml(calc_input_file, lj_interactive_calc_input)
    write_yaml(sim_input_file, singlepoint_sim_input)

    asim_run([calc_input_file, sim_input_file])

    # Check job completion
    output = read_yaml('test_output.yaml')
    assert output.get('status', False) == 'complete'

def test_asim_run_do_nothing(
    tmp_path,
    lj_interactive_calc_input,
    do_nothing_sim_input,
):
    ''' Test running an arbitrary script defined by user '''
    # Create input files in temporary directory that will 
    # be automatically deleted using the fixtures in conftest.py
    os.chdir (tmp_path)
    calc_input_file = 'calc_input.yaml'
    sim_input_file = 'sim_input.yaml'
    write_yaml(calc_input_file, lj_interactive_calc_input)

    sim_input = do_nothing_sim_input.copy()
    sim_input['script'] = SCRIPT_DIR
    write_yaml(sim_input_file, sim_input)

    asim_run([calc_input_file, sim_input_file])

    # Check job completion
    output = read_yaml('output.yaml')
    assert output.get('status', False) == 'complete'
