'''
Test Job class
'''
#pylint: disable=missing-function-docstring
#pylint: disable=redefined-outer-name

import pytest
from asimtools.job import UnitJob

def create_unitjob(sim_input, env_input, workdir, calc_input=None):
    """Helper for making a generic UnitJob object"""
    env_id = list(env_input.keys())[0]
    sim_input['env_id'] = env_id
    if calc_input is not None:
        calc_id = list(calc_input.keys())[0]
        sim_input['calc_id'] = calc_id
    sim_input['workdir'] = workdir
    unitjob = UnitJob(
        sim_input,
        env_input,
        calc_input
    )
    return unitjob

@pytest.mark.parametrize("env_input",[
    "inline_env_input", "batch_env_input", "salloc_env_input"
])
@pytest.mark.parametrize("sim_input",["do_nothing_sim_input"])
def test_gen_input_files(env_input, sim_input, tmp_path, request):
    env_input = request.getfixturevalue(env_input)
    sim_input = request.getfixturevalue(sim_input)
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(sim_input, env_input, wdir)
    unitjob.gen_input_files()

    assert wdir.exists()
    assert (wdir / 'sim_input.yaml').exists()
    assert (wdir / 'calc_input.yaml').exists()
    assert (wdir / 'env_input.yaml').exists()
    assert (wdir / 'output.yaml').exists()

@pytest.mark.parametrize("env_input",[
    "inline_env_input", "batch_env_input", "salloc_env_input"
])
@pytest.mark.parametrize("sim_input",["do_nothing_sim_input"])
def test_update_and_read_output(env_input, sim_input, tmp_path, request):
    env_input = request.getfixturevalue(env_input)
    sim_input = request.getfixturevalue(sim_input)
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(sim_input, env_input, wdir)
    unitjob.gen_input_files()

    output_update = {'test_val': 1}
    unitjob.gen_input_files()
    unitjob.update_output(output_update)
    assert unitjob.get_output().get('test_val', False) == 1

@pytest.mark.parametrize("env_input",[
    "inline_env_input", "batch_env_input", "salloc_env_input"
])
@pytest.mark.parametrize("sim_input",["do_nothing_sim_input"])
def test_start_fail_discard_complete(env_input, sim_input, tmp_path, request):
    env_input = request.getfixturevalue(env_input)
    sim_input = request.getfixturevalue(sim_input)
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(sim_input, env_input, wdir)
    unitjob.gen_input_files()

    unitjob.start()
    assert unitjob.get_output().get('start_time', False)
    assert unitjob.get_status() == (False, 'started')
    unitjob.fail()
    assert unitjob.get_status() == (False, 'failed')
    unitjob.discard()
    assert unitjob.get_status() == (False, 'discard')
    unitjob.complete()
    assert unitjob.get_output().get('end_time', False)
    assert unitjob.get_status() == (True, 'complete')

@pytest.mark.parametrize("env_input",["inline_env_input"])
@pytest.mark.parametrize("sim_input",["do_nothing_sim_input"])
@pytest.mark.parametrize("submit, overwrite, expected",[
    (True, True, (True, 'complete')),
    (False, True, (False, 'clean')),
    (True, False, (True, 'complete')),
    (False, False, (False, 'clean')),
])
def test_submit(
    env_input, sim_input, submit, overwrite, expected, tmp_path, request
):
    env_input = request.getfixturevalue(env_input)
    sim_input = request.getfixturevalue(sim_input)
    wdir = tmp_path / 'wdir'

    # Prepare but don't submit
    sim_input['submit'] = submit
    sim_input['overwrite'] = overwrite
    unitjob = create_unitjob(sim_input, env_input, wdir)
    unitjob.submit()
    assert unitjob.get_status() == expected

    unitjob = create_unitjob(sim_input, env_input, wdir)
    if submit:
        old_start_time = unitjob.get_output()['start_time']
        unitjob.submit()    
        new_start_time = unitjob.get_output()['start_time']
        if overwrite is True:
            assert new_start_time != old_start_time
        else:
            assert new_start_time == old_start_time
