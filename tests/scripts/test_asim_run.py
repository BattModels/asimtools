"""
Tests for running asimmodules using asim_run.py
"""
import os
from pathlib import Path
import pytest
from asimtools.utils import write_yaml, read_yaml
from asimtools.job import create_unitjob
from asimtools.scripts.asim_run import main as asim_run
from asimtools.scripts.asim_run import parse_command_line

@pytest.mark.parametrize("calc_input",["lj_argon_calc_input"])
@pytest.mark.parametrize("sim_input",["do_nothing_sim_input"])
def test_parse_command_line(sim_input, calc_input, tmp_path):
    ''' Tests commandline argument parsing '''
    calc_input_file = str(tmp_path / 'calc_input.yaml')
    sim_input_file = str(tmp_path / 'sim_input.yaml')
    write_yaml(calc_input_file, calc_input)
    write_yaml(sim_input_file, sim_input)

    sim_config_parse, calc_config_parse = parse_command_line(
        [sim_input_file, '-c', calc_input_file]
    )

    assert Path(calc_input_file).resolve() == Path(calc_config_parse).resolve()
    assert sim_input == sim_config_parse

@pytest.mark.parametrize("inputs",[
    ['do_nothing_sim_input', None],
    ['do_nothing_sim_input', 'lj_argon_calc_input'],
])
def test_asim_run(inputs, tmp_path, request):
    ''' 
    Test asim_run ability to launch and complete a asimmodule using LJ
    '''
    os.chdir(tmp_path)
    sim_input = request.getfixturevalue(inputs[0])
    sim_input_file = 'sim_input.yaml'
    write_yaml(sim_input_file, sim_input)
    args = [sim_input_file]

    if inputs[1] is not None:
        calc_input_file = 'calc_input.yaml'
        calc_input = request.getfixturevalue(inputs[1])
        write_yaml(calc_input_file, calc_input)
        args += ['-c', calc_input_file]

    asim_run(args)

    # Check job completion
    output = read_yaml('output.yaml')
    assert output.get('status', False) == 'complete'

@pytest.mark.parametrize("inputs",[
    ['do_nothing_sim_input', None],
])
def test_asim_run_commands(inputs, tmp_path, request):
    ''' 
    Test asim_run ability to launch and complete a asimmodule using LJ and running
    commands before and after running the function
    '''
    os.chdir(tmp_path)
    sim_input = request.getfixturevalue(inputs[0])
    sim_input['precommands'] = ['touch precommand_test']
    sim_input['postcommands'] = ['touch postcommand_test']
    sim_input_file = 'sim_input.yaml'
    write_yaml(sim_input_file, sim_input)
    args = [sim_input_file]

    if inputs[1] is not None:
        calc_input_file = 'calc_input.yaml'
        calc_input = request.getfixturevalue(inputs[1])
        write_yaml(calc_input_file, calc_input)
        args += ['-c', calc_input_file]

    asim_run(args)

    assert Path('precommand_test').exists()
    assert Path('postcommand_test').exists()
