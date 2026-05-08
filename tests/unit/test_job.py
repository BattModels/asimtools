'''
Test Job class
'''
#pylint: disable=missing-function-docstring
#pylint: disable=redefined-outer-name

import logging
import os
from pathlib import Path
import pytest
from asimtools.job import (
    UnitJob,
    DistributedJob,
    ChainedJob,
    check_job_tree_complete,
    load_job_from_directory,
    get_subjobs,
    load_job_tree,
)
from asimtools.utils import read_yaml, write_yaml

def create_unitjob(sim_input, env_input, workdir, calc_input=None, status=None):
    """Helper for making a generic UnitJob object"""
    env_id = list(env_input.keys())[0]
    sim_input['env_id'] = env_id
    if calc_input is not None:
        calc_id = list(calc_input.keys())[0]
        sim_input.setdefault('args', {})['calculator'] = {'calc_id': calc_id}
    sim_input['workdir'] = workdir
    unitjob = UnitJob(
        sim_input,
        env_input,
        calc_input
    )
    if status is not None:
        unitjob.set_status(status)
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
        'args': {'calculator': {'calc_id': 'test_calc_id'}},
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

class FloatingJob(UnitJob):
    '''A job that can be in any status without being written to disk'''
    def __init__(self, sim_input, env_input, calc_input=None, status='clean'):
        super().__init__(sim_input, env_input, calc_input)
        self.status = status

    def get_status(self):
        completed = False
        if self.status == 'complete':
            completed = True
        return completed, self.status

complete_job = FloatingJob(
    {'asimmodule': 'do_nothing'},
    {'inline_env': {'mode': {'use_slurm': False}}},
    status='complete',  
)

discard_job = FloatingJob(
    {'asimmodule': 'do_nothing'},
    {'inline_env': {'mode': {'use_slurm': False}}},
    status='discard',
)

failed_job = FloatingJob(
    {'asimmodule': 'do_nothing'},
    {'inline_env': {'mode': {'use_slurm': False}}},
    status='failed',
)

clean_job = FloatingJob(
    {'asimmodule': 'do_nothing'},
    {'inline_env': {'mode': {'use_slurm': False}}},
    status='clean',
)

@pytest.mark.parametrize("test_input, expected",[
    ({
        'workdir_name': 'dummy',
        'job': complete_job,
        'subjobs': None,
    }, (True, 'complete')),
    ({
        'workdir_name': 'dummy',
        'job': failed_job,
        'subjobs': None,
    }, (False, 'failed')),
    ({
        'workdir_name': 'dummy',
        'job': clean_job,
        'subjobs': None,
    }, (False, 'clean')),
    ({
        'workdir_name': 'dummy',
        'job': discard_job,
        'subjobs': None,
    }, (True, 'discard')),
    ({
        'workdir_name': 'dummy',
        'job': complete_job,
        'subjobs': {
            'id-0000': {
                'workdir_name': 'dummy',
                'job': complete_job,
                'subjobs': None,
            },
            'id-0001': {
                'workdir_name': 'dummy',
                'job': complete_job,
                'subjobs': None,
            },
        }
    }, (True, 'complete')),
    ({
        'workdir_name': 'dummy',
        'job': complete_job,
        'subjobs': {
            'id-0000': {
                'workdir_name': 'dummy',
                'job': discard_job,
                'subjobs': None,
            },
            'id-0001': {
                'workdir_name': 'dummy',
                'job': complete_job,
                'subjobs': None,
            },
        }
    }, (True, 'complete')),
    ({
        'workdir_name': 'dummy',
        'job': complete_job,
        'subjobs': {
            'id-0000': {
                'workdir_name': 'dummy',
                'job': complete_job,
                'subjobs': None,
            },
            'id-0001': {
                'workdir_name': 'dummy',
                'job': failed_job,
                'subjobs': None,
            },
        }
    }, (False, 'failed')),
    ({
        'workdir_name': 'dummy',
        'job': complete_job,
        'subjobs': {
            'id-0000': {
                'workdir_name': 'dummy',
                'job': complete_job,
                'subjobs': None,
            },
            'id-0001': {
                'workdir_name': 'dummy',
                'job': complete_job,
                'subjobs': {
                    'id-0000': {
                        'workdir_name': 'dummy',
                        'job': discard_job,
                        'subjobs': None,
                    },
                    'id-0001': {
                        'workdir_name': 'dummy',
                        'job': clean_job,
                        'subjobs': None,
                    },
                },
            },
        }
    }, (False, 'clean')),
    ({
        'workdir_name': 'dummy',
        'job': complete_job,
        'subjobs': {
            'id-0000': {
                'workdir_name': 'dummy',
                'job': complete_job,
                'subjobs': None,
            },
            'id-0001': {
                'workdir_name': 'dummy',
                'job': complete_job,
                'subjobs': {
                    'id-0000': {
                        'workdir_name': 'dummy',
                        'job': complete_job,
                        'subjobs': None,
                    },
                    'id-0001': {
                        'workdir_name': 'dummy',
                        'job': complete_job,
                        'subjobs': None,
                    },
                },
            },
        }
    }, (True, 'complete')),
])
def test_check_job_tree_complete(tmp_path, test_input, expected):
    ''' Test check_job_tree_complete '''

    assert check_job_tree_complete(test_input) == expected
    assert check_job_tree_complete(test_input, skip_failed=True)[0] is True


# ---------------------------------------------------------------------------
# Job getter / updater methods
# ---------------------------------------------------------------------------

def test_get_sim_input_returns_internal(inline_env_input, do_nothing_sim_input, tmp_path):
    unitjob = create_unitjob(do_nothing_sim_input, inline_env_input, tmp_path / 'wdir')
    assert unitjob.get_sim_input() is unitjob.sim_input


def test_get_calc_input_returns_internal(inline_env_input, do_nothing_sim_input, lj_argon_calc_input, tmp_path):
    unitjob = create_unitjob(
        do_nothing_sim_input, inline_env_input, tmp_path / 'wdir', calc_input=lj_argon_calc_input
    )
    assert unitjob.get_calc_input() is unitjob.calc_input


def test_get_env_input_returns_internal(inline_env_input, do_nothing_sim_input, tmp_path):
    unitjob = create_unitjob(do_nothing_sim_input, inline_env_input, tmp_path / 'wdir')
    assert unitjob.get_env_input() is unitjob.env_input


def test_update_sim_input(inline_env_input, do_nothing_sim_input, tmp_path):
    unitjob = create_unitjob(do_nothing_sim_input, inline_env_input, tmp_path / 'wdir')
    unitjob.update_sim_input({'new_key': 'new_value'})
    assert unitjob.get_sim_input()['new_key'] == 'new_value'


def test_add_output_files(inline_env_input, do_nothing_sim_input, tmp_path):
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(do_nothing_sim_input, inline_env_input, wdir)
    unitjob.gen_input_files()
    unitjob.add_output_files({'result': 'out.xyz'})
    assert unitjob.get_output().get('files', {}).get('result') == 'out.xyz'


def test_get_logger(inline_env_input, do_nothing_sim_input, tmp_path):
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(do_nothing_sim_input, inline_env_input, wdir)
    unitjob.gen_input_files()
    logger = unitjob.get_logger()
    assert isinstance(logger, logging.Logger)


# ---------------------------------------------------------------------------
# UnitJob.gen_run_command
# ---------------------------------------------------------------------------

def test_gen_run_command_basic(inline_env_input, do_nothing_sim_input, tmp_path):
    unitjob = create_unitjob(do_nothing_sim_input, inline_env_input, tmp_path / 'wdir')
    cmd = unitjob.gen_run_command()
    assert 'asim-run' in cmd
    assert 'sim_input.yaml' in cmd


def test_gen_run_command_prefix_suffix(tmp_path):
    env_input = {
        'inline': {
            'mode': {'use_slurm': False, 'interactive': True},
        }
    }
    calc_input = {
        'lj': {
            'name': 'LennardJones',
            'module': 'ase.calculators.lj',
            'run_prefix': 'mpirun -n 4',
            'run_suffix': '--some-flag',
            'args': {},
        }
    }
    sim_input = {
        'asimmodule': 'do_nothing',
        'env_id': 'inline',
        'workdir': str(tmp_path / 'wdir'),
        'args': {'calculator': {'calc_id': 'lj'}},
    }
    unitjob = UnitJob(sim_input, env_input=env_input, calc_input=calc_input)
    cmd = unitjob.gen_run_command()
    assert cmd.index('mpirun -n 4') < cmd.index('asim-run')
    assert cmd.index('asim-run') < cmd.index('--some-flag')


# ---------------------------------------------------------------------------
# DistributedJob
# ---------------------------------------------------------------------------

def _make_distributed_sim_input(env_id, n=3):
    subsim = {'asimmodule': 'do_nothing', 'env_id': env_id}
    return {f'job{i}': dict(subsim) for i in range(n)}


def test_distributed_job_init_inline(inline_env_input):
    sim_input = _make_distributed_sim_input('inline', n=3)
    djob = DistributedJob(sim_input, env_input=inline_env_input)
    assert len(djob.unitjobs) == 3
    assert djob.use_slurm is False
    # Workdir names must be prefixed with 'id-'
    for uj in djob.unitjobs:
        assert uj.workdir.name.startswith('id-')


def test_distributed_job_init_slurm(batch_env_input):
    sim_input = _make_distributed_sim_input('batch', n=2)
    djob = DistributedJob(sim_input, env_input=batch_env_input)
    assert djob.use_slurm is True


def test_distributed_job_init_too_many_jobs(inline_env_input):
    sim_input = {f'job{i}': {'asimmodule': 'do_nothing', 'env_id': 'inline'} for i in range(1000)}
    with pytest.raises(AssertionError):
        DistributedJob(sim_input, env_input=inline_env_input)


def test_distributed_job_gen_input_files(inline_env_input, tmp_path):
    original_dir = Path('.').resolve()
    os.chdir(tmp_path)
    try:
        sim_input = _make_distributed_sim_input('inline', n=2)
        djob = DistributedJob(sim_input, env_input=inline_env_input)
        djob.gen_input_files()
        for uj in djob.unitjobs:
            assert (uj.workdir / 'sim_input.yaml').exists()
    finally:
        os.chdir(original_dir)


def test_distributed_job_submit_inline(inline_env_input, tmp_path):
    original_dir = Path('.').resolve()
    os.chdir(tmp_path)
    try:
        sim_input = _make_distributed_sim_input('inline', n=2)
        djob = DistributedJob(sim_input, env_input=inline_env_input)
        djob.submit()
        for uj in djob.unitjobs:
            assert uj.get_status()[1] == 'complete'
    finally:
        os.chdir(original_dir)


# ---------------------------------------------------------------------------
# ChainedJob
# ---------------------------------------------------------------------------

def _make_chained_sim_input(env_id, n=2):
    return {f'step-{i}': {'asimmodule': 'do_nothing', 'env_id': env_id} for i in range(n)}


def test_chained_job_init(inline_env_input):
    sim_input = _make_chained_sim_input('inline', n=3)
    cjob = ChainedJob(sim_input, env_input=inline_env_input)
    assert len(cjob.unitjobs) == 3
    for i, uj in enumerate(cjob.unitjobs):
        assert str(uj.workdir) == f'step-{i}'


def test_chained_job_init_bad_keys(inline_env_input):
    sim_input = {'bad_key': {'asimmodule': 'do_nothing', 'env_id': 'inline'}}
    with pytest.raises(AssertionError):
        ChainedJob(sim_input, env_input=inline_env_input)


def test_chained_job_get_last_output(tmp_path):
    sim_input = _make_chained_sim_input('inline', n=2)
    env_input = {'inline': {'mode': {'use_slurm': False, 'interactive': True}}}
    chain_dir = tmp_path / 'chain'
    chain_dir.mkdir()
    cjob = ChainedJob(sim_input, env_input=env_input)
    cjob.set_workdir(chain_dir)
    cjob.unitjobs[-1].set_workdir(chain_dir / 'step-1')
    cjob.unitjobs[-1].gen_input_files()
    cjob.unitjobs[-1].update_output({'my_result': 42})
    assert cjob.get_last_output()['my_result'] == 42


# ---------------------------------------------------------------------------
# load_job_from_directory, get_subjobs, load_job_tree
# ---------------------------------------------------------------------------

def test_load_job_from_directory(inline_env_input, do_nothing_sim_input, tmp_path):
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(do_nothing_sim_input, inline_env_input, wdir)
    unitjob.gen_input_files(write_env_input=True)
    job = load_job_from_directory(wdir)
    assert job.sim_input['asimmodule'] == do_nothing_sim_input['asimmodule']
    assert job.workdir == wdir


def test_load_job_from_directory_missing(tmp_path):
    with pytest.raises(AssertionError):
        load_job_from_directory(tmp_path / 'nonexistent')


def test_get_subjobs(inline_env_input, do_nothing_sim_input, tmp_path):
    root = tmp_path / 'root'
    root.mkdir()
    write_yaml(root / 'sim_input.yaml', do_nothing_sim_input)

    # Two subdirs with sim_input.yaml, one without
    for name in ('sub_a', 'sub_b', 'empty_sub'):
        (root / name).mkdir()
    write_yaml(root / 'sub_a' / 'sim_input.yaml', do_nothing_sim_input)
    write_yaml(root / 'sub_b' / 'sim_input.yaml', do_nothing_sim_input)

    subjobs = get_subjobs(root)
    assert len(subjobs) == 2
    assert all(p.name in ('sub_a', 'sub_b') for p in subjobs)
    # Must be sorted
    assert subjobs == sorted(subjobs)


def test_load_job_tree_flat(inline_env_input, do_nothing_sim_input, tmp_path):
    wdir = tmp_path / 'wdir'
    unitjob = create_unitjob(do_nothing_sim_input, inline_env_input, wdir)
    unitjob.gen_input_files(write_env_input=True)
    tree = load_job_tree(wdir)
    assert tree['workdir_name'] == wdir.name
    assert tree['subjobs'] is None


def test_load_job_tree_nested(inline_env_input, do_nothing_sim_input, tmp_path):
    root = tmp_path / 'root'
    root.mkdir()
    write_yaml(root / 'sim_input.yaml', do_nothing_sim_input)
    for name in ('sub_a', 'sub_b'):
        subdir = root / name
        subdir.mkdir()
        write_yaml(subdir / 'sim_input.yaml', do_nothing_sim_input)

    tree = load_job_tree(root)
    assert tree['subjobs'] is not None
    assert set(tree['subjobs'].keys()) == {'sub_a', 'sub_b'}
    assert tree['subjobs']['sub_a']['subjobs'] is None
