"""
Tests for running scripts using asim_run.py
"""
import os
import glob
from pathlib import Path
import pytest
from asimtools.utils import write_yaml
from asimtools.scripts.asim_execute import main as asim_execute

@pytest.mark.parametrize("inputs",[
    ['do_nothing_distributed_sim_input', None],
])
def test_distributed(inputs, tmp_path, request):
    ''' 
    Test distributed ability to launch and complete distributed jobs
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

    asim_execute(args)
    print(str(Path(tmp_path) / 'id-00')+'*')
    print(glob.glob(str(Path(tmp_path))+'/*'))

    assert len(glob.glob(str(Path(tmp_path) / 'id-00')+'*')) == 1
