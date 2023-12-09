'''
Test Job class
'''
#pylint: disable=missing-function-docstring
#pylint: disable=redefined-outer-name

from pathlib import Path
import pytest
from asimtools.job import UnitJob
from asimtools.utils import read_yaml

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
    "inline_env_input",
    "batch_env_input",
    "salloc_env_input",
])
@pytest.mark.parametrize("sim_input",[
    "do_nothing_sim_input", "singlepoint_argon_sim_input"
])
def test_gen_input_files(env_input, sim_input, tmp_path, request):
    env_input = request.getfixturevalue(env_input)
    sim_input = request.getfixturevalue(sim_input)
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(sim_input, env_input, wdir)
    unitjob.gen_input_files(write_env_input=True, write_calc_input=True)

    assert wdir.exists()
    assert (wdir / 'sim_input.yaml').exists()
    assert (wdir / 'calc_input.yaml').exists()
    assert (wdir / 'env_input.yaml').exists()
    assert (wdir / 'output.yaml').exists()
    assert Path(read_yaml(wdir / 'sim_input.yaml')['workdir']) == Path('./')

    # Test that we always write the input atomic structures
    if "image" in sim_input.get('args', {}).get('image', []):
        assert (wdir / 'input_image.xyz').exists()
    if "images" in sim_input.get('args', {}).get('images', []):
        assert (wdir / 'input_images.xyz').exists()

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

def test_slurm_script(tmp_path):  
    wdir = tmp_path / 'results'
    sim_input = {
        'script': 'do_nothing',
        'env_id': 'test_batch',
        'workdir': wdir,
        'args': {'calc_id': 'test_calc_id'},
    }
    env_input = {
        'test_batch': {
            'mode': {
                'use_slurm': True,
                'interactive': False,
            },
            'slurm': {
                'flags': ['-N 1 #test flag'],
                'precommands': ['test_env_precommand1', 'test_env_precommand2'],
                'postcommands': ['test_env_postcommand1', 'test_env_postcommand2'],
            },
        },
    }

    calc_input = {
        'test_calc_id': {
            'precommands': ['test_calc_precommand1', 'test_calc_precommand2'],
            'postcommands': ['test_calc_postcommand1', 'test_calc_postcommand2'],
            'run_prefix': 'touch test_run_prefix; ',
            'run_suffix': '; touch test_run_suffix',
        },
    }

    unitjob = create_unitjob(sim_input, env_input, wdir, calc_input=calc_input)
    unitjob.gen_input_files()

    with open(wdir / 'job.sh', 'r', encoding='utf-8') as f:
        lines = f.readlines()
    jscript_txt = ' '.join(lines)

    assert '#SBATCH -N 1' in jscript_txt
    assert 'test_env_precommand1' in jscript_txt
    assert 'test_env_precommand2' in jscript_txt
    assert 'test_env_postcommand1' in jscript_txt
    assert 'test_env_postcommand2' in jscript_txt
    assert 'test_calc_precommand1' in jscript_txt
    assert 'test_calc_precommand2' in jscript_txt
    assert 'test_calc_postcommand1' in jscript_txt
    assert 'test_calc_postcommand2' in jscript_txt
    assert 'test_run_prefix' in jscript_txt
    assert 'test_run_suffix' in jscript_txt

    for line in lines:
        if 'test_run_prefix' in line:
            assert 'asim-run' in line
            assert line.index('test_run_prefix') < line.index('asim-run')
        if 'test_run_suffix' in line:
            assert 'asim-run' in line
            assert line.index('test_run_suffix') > line.index('asim-run')
