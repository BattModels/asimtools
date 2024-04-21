"""
Tests for running asimmodules using asim_run.py
"""
from glob import glob
import pytest
from asimtools.job import create_unitjob
from asimtools.job import load_job_from_directory

@pytest.mark.parametrize("calc_input",["lj_argon_calc_input"])
@pytest.mark.parametrize("env_input",["inline_env_input"])
@pytest.mark.parametrize("sim_input",["lj_distributed_sim_input"])
def test_distributed(env_input, calc_input, sim_input, tmp_path, request):
    env_input = request.getfixturevalue(env_input)
    calc_input = request.getfixturevalue(calc_input)
    sim_input = request.getfixturevalue(sim_input)
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(sim_input, env_input, wdir, calc_input=calc_input)
    unitjob.submit()

    assert load_job_from_directory(wdir).get_status()[1] == 'complete'
    dirs = glob(str(wdir / 'id*'))
    assert len(dirs) == len(sim_input['args']['subsim_inputs'])

    for d in dirs:
        assert str(d).rsplit('/', maxsplit=1)[-1].startswith('id-')

        uj = load_job_from_directory(d)
        assert uj.get_status()[1] == 'complete'

        assert uj.get_sim_input()['workdir'] == './'
