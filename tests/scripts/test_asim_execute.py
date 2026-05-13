"""
Tests for running asimmodules using asim_execute.py
"""
import pytest
from asimtools.job import create_unitjob


@pytest.mark.parametrize("env_input", ["inline_env_input"])
@pytest.mark.parametrize("sim_input", ["do_nothing_sim_input"])
def test_asim_execute_dry_run(sim_input, env_input, tmp_path, request):
    """dry_run=True writes input files but does not submit the job"""
    sim_input = request.getfixturevalue(sim_input)
    env_input = request.getfixturevalue(env_input)
    sim_input['dry_run'] = True
    wdir = tmp_path / 'wdir'

    unitjob = create_unitjob(sim_input, env_input, wdir)
    unitjob.submit()

    assert (wdir / 'sim_input.yaml').exists()
    assert unitjob.get_status()[1] == 'clean'
