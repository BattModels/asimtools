"""
Tests for running asimmodules using asim_run.py
"""
from glob import glob
import os
import pytest
from natsort import natsorted
from asimtools.job import (
    create_unitjob,
    DistributedJob,
    load_job_from_directory,
)

def create_distjob(sim_input, env_input, workdir, calc_input=None):
    """Helper for making a generic DistributedJob object, mostly for testing"""
    env_id = list(env_input.keys())[0]
    sim_input['env_id'] = env_id
    if calc_input is not None:
        calc_id = list(calc_input.keys())[0]
        sim_input['calc_id'] = calc_id
    sim_input['workdir'] = workdir
    distjob = DistributedJob(
        sim_input['args']['subsim_inputs'],
        env_input,
        calc_input,
    )
    return distjob

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
    dirs = natsorted(glob(str(wdir / 'id*')))
    assert len(dirs) == len(sim_input['args']['subsim_inputs'])

    for d in dirs:
        assert str(d).rsplit('/', maxsplit=1)[-1].startswith('id-')

        uj = load_job_from_directory(d)
        assert uj.get_status()[1] == 'complete'

        assert uj.get_sim_input()['workdir'] == './'

    assert unitjob.get_status(descend=True) == (True, 'complete')

@pytest.mark.parametrize("calc_input",["lj_argon_calc_input"])
@pytest.mark.parametrize("env_input",["batch_env_input"])
@pytest.mark.parametrize("sim_input",[
    "lj_distributed_batch_sim_input",
    "lj_distributed_group_batch_sim_input",
])
def test_batch_distributed(env_input, calc_input, sim_input, tmp_path, request):
    env_input = request.getfixturevalue(env_input)
    calc_input = request.getfixturevalue(calc_input)
    sim_input = request.getfixturevalue(sim_input)
    group_size = sim_input['args'].get('group_size', 1)
    wdir = tmp_path
    os.chdir(wdir)
    distjob = create_distjob(sim_input, env_input, wdir, calc_input=calc_input)
    group_size=sim_input['args'].get('group_size', 1)
    array_max=sim_input['args'].get('array_max', None)
    skip_failed=sim_input['args'].get('skip_failed', False)
    distjob.submit(
        debug=True,
        group_size=group_size,
        array_max=array_max,
        skip_failed=skip_failed,
    )
    dirs = natsorted(glob(str(wdir / 'id*')))
    assert len(dirs) == len(sim_input['args']['subsim_inputs'])

    # In debug mode for distjob, only the first group will be run
    statuses = ['complete'] * group_size + ['clean'] * (len(dirs) - 1)
    print("XXXXXXX", statuses, dirs)
    for d_ind, d in enumerate(dirs):
        assert str(d).rsplit('/', maxsplit=1)[-1].startswith('id-')

        uj = load_job_from_directory(d)
        print('job_info:', uj.workdir, uj.get_status())
        assert uj.get_status()[1] == statuses[d_ind]

    # assert distjob.get_status(descend=False) == (True, 'complete')
