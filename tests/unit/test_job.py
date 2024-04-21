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
    "batch_dict_env_input",
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
        assert (wdir / 'image_input.xyz').exists()
    if "images" in sim_input.get('args', {}).get('images', []):
        assert (wdir / 'images_input.xyz').exists()

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
@pytest.mark.parametrize("calc_input",["lj_argon_calc_input"])
@pytest.mark.parametrize("sim_input",["do_nothing_sim_input"])
def test_update_calc_env_input(env_input, sim_input, calc_input, tmp_path, request):
    env_input = request.getfixturevalue(env_input)
    sim_input = request.getfixturevalue(sim_input)
    calc_input = request.getfixturevalue(calc_input)
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(sim_input, env_input, wdir)
    env_id = unitjob.env_id

    new_env_input = env_input.copy()
    new_env_input[env_id]['precommands'] = ['precommand']
    unitjob.update_env_input(new_env_input)

    assert 'precommand' in unitjob.env_input[env_id]['precommands']
    assert 'precommand' in unitjob.env['precommands']

    new_calc_input = calc_input.copy()
    calc_id = list(calc_input.keys())[0]
    new_calc_input[calc_id]['precommands'] = ['precommand']
    unitjob.update_calc_input(new_calc_input)
    assert 'precommand' in unitjob.calc_input[calc_id]['precommands']

@pytest.mark.parametrize("env_input",[
    "inline_env_input",
    "batch_env_input",
    "batch_dict_env_input",
    "salloc_env_input",
])
@pytest.mark.parametrize("sim_input",[
    "do_nothing_sim_input", "singlepoint_argon_sim_input"
])
def test_set_workdir(env_input, sim_input, tmp_path, request):
    env_input = request.getfixturevalue(env_input)
    sim_input = request.getfixturevalue(sim_input)
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(sim_input, env_input, wdir)

    new_workdir = (tmp_path / 'new_wdir').resolve()
    unitjob.set_workdir(new_workdir)
    assert unitjob.get_workdir().resolve() == new_workdir
    assert unitjob.workdir.resolve() == new_workdir
    assert Path(unitjob.sim_input['workdir']).resolve() == new_workdir
    assert isinstance(unitjob.sim_input['workdir'], str)

    str_workdir = str((tmp_path / 'new_wdir').resolve())
    unitjob.set_workdir(str_workdir)
    assert unitjob.get_workdir().resolve() == new_workdir
    assert unitjob.workdir.resolve() == new_workdir
    assert Path(unitjob.sim_input['workdir']).resolve() == new_workdir

@pytest.mark.parametrize("env_input",[
    "inline_env_input", "batch_env_input", "batch_dict_env_input", "salloc_env_input"
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

@pytest.mark.parametrize("env_input",[
    "inline_env_input", "batch_env_input", "batch_dict_env_input", "salloc_env_input"
])
@pytest.mark.parametrize("sim_input",["do_nothing_sim_input"])
def test_go_to_leave_workdir(env_input, sim_input, tmp_path, request):
    env_input = request.getfixturevalue(env_input)
    sim_input = request.getfixturevalue(sim_input)
    cur_dir = Path('.').resolve()
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(sim_input, env_input, wdir)
    unitjob.gen_input_files()
    unitjob.go_to_workdir()
    assert Path('.').resolve() == wdir.resolve()
    unitjob.leave_workdir()
    assert Path('.').resolve() == cur_dir

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

@pytest.mark.parametrize("flags",(
    ['-N 1 #test flag', '-n 4', '--mem-per-cpu=4G'],
    {'-N': 1 , '-n': '4', '--mem-per-cpu': '4G'},
))
def test_slurm_asimmodule(flags, tmp_path):  
    wdir = tmp_path / 'results'
    jobname = 'test_job'
    sim_input = {
        'asimmodule': 'do_nothing',
        'env_id': 'test_batch',
        'workdir': wdir,
        'job_name': jobname,
        'args': {'calc_id': 'test_calc_id'},
    }

    env_input = {
        'test_batch': {
            'mode': {
                'use_slurm': True,
                'interactive': False,
            },
            'slurm': {
                'flags': flags,
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
    jasimmodule_txt = ' '.join(lines)

    assert '#SBATCH -N 1' in jasimmodule_txt
    assert '#SBATCH -n 4' in jasimmodule_txt
    assert f'#SBATCH -J {jobname}' in jasimmodule_txt
    assert '#SBATCH --mem-per-cpu=4' in jasimmodule_txt
    assert 'test_env_precommand1' in jasimmodule_txt
    assert 'test_env_precommand2' in jasimmodule_txt
    assert 'test_env_postcommand1' in jasimmodule_txt
    assert 'test_env_postcommand2' in jasimmodule_txt
    assert 'test_calc_precommand1' in jasimmodule_txt
    assert 'test_calc_precommand2' in jasimmodule_txt
    assert 'test_calc_postcommand1' in jasimmodule_txt
    assert 'test_calc_postcommand2' in jasimmodule_txt
    assert 'test_run_prefix' in jasimmodule_txt
    assert 'test_run_suffix' in jasimmodule_txt

    for line in lines:
        if 'test_run_prefix' in line:
            assert 'asim-run' in line
            assert line.index('test_run_prefix') < line.index('asim-run')
        if 'test_run_suffix' in line:
            assert 'asim-run' in line
            assert line.index('test_run_suffix') > line.index('asim-run')
