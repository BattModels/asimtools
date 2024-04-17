"""
Tests for running asimmodules using asim_run.py
"""
from pathlib import Path
from glob import glob
import pytest
from asimtools.utils import write_yaml
from asimtools.job import create_unitjob
from asimtools.job import load_job_from_directory
from asimtools.scripts.asim_check import (
    get_subjobs,
    load_job_tree,
)
from asimtools.scripts.asim_check import parse_command_line

@pytest.mark.parametrize("sim_input",["do_nothing_sim_input"])
def test_parse_command_line(sim_input, tmp_path):
    ''' Tests commandline argument parsing '''
    sim_input_file = str(tmp_path / 'sim_input.yaml')
    write_yaml(sim_input_file, sim_input)

    sim_config_parse, rootdir, max_level = parse_command_line(
        [sim_input_file, '-m 2']
    )

    assert sim_input == sim_config_parse
    assert Path(rootdir) == Path(tmp_path)
    assert max_level == 2

@pytest.mark.parametrize("calc_input",["lj_argon_calc_input"])
@pytest.mark.parametrize("env_input",["inline_env_input"])
@pytest.mark.parametrize("sim_input",["lj_distributed_sim_input"])
def test_get_subjobs(env_input, calc_input, sim_input, tmp_path, request):
    env_input = request.getfixturevalue(env_input)
    calc_input = request.getfixturevalue(calc_input)
    sim_input = request.getfixturevalue(sim_input)
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(sim_input, env_input, wdir, calc_input=calc_input)
    unitjob.submit()

    dirs = glob(str(wdir / 'id*'))
    assert len(dirs) == len(sim_input['args']['subsim_inputs'])

    ac_dirs = get_subjobs(wdir)
    assert len(ac_dirs) == len(sim_input['args']['subsim_inputs'])
    for d in ac_dirs:
        assert load_job_from_directory(d).get_status()[1] == 'complete'

@pytest.mark.parametrize("calc_input",["lj_argon_calc_input"])
@pytest.mark.parametrize("env_input",["inline_env_input"])
@pytest.mark.parametrize("sim_input",["lj_distributed_sim_input"])
def test_load_job_tree(env_input, calc_input, sim_input, tmp_path, request):
    env_input = request.getfixturevalue(env_input)
    calc_input = request.getfixturevalue(calc_input)
    sim_input = request.getfixturevalue(sim_input)
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(sim_input, env_input, wdir, calc_input=calc_input)
    unitjob.submit()

    jd = load_job_tree(wdir)
    print(jd)
    assert jd['workdir_name'] == 'wdir'
    assert len(jd['subjobs']) == len(sim_input['args']['subsim_inputs']) 
    for sj in jd['subjobs']:
        assert jd['subjobs'][sj]['subjobs'] is None
        assert jd['subjobs'][sj].get('workdir_name', False)
