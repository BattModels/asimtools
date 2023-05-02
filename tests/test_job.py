'''
Test Job class
'''
#pylint: disable=missing-function-docstring
#pylint: disable=redefined-outer-name

import builtins
import os
import io
import pytest
from asimtools.job import UnitJob

# def patch_open(open_func, files):
#     '''
#     For patching open so that we clean up any files we write to disk 
#     after tests. See:
#     https://stackoverflow.com/questions/51737334/pytest-deleting-
#     files-created-by-the-tested-function
#     '''
#     def open_patched(path, mode='r', buffering=-1, encoding=None, 
#                     errors=None, newline=None, closefd=True,
#                     opener=None):
#         if 'w' in mode and not os.path.isfile(path):
#             files.append(path)
#         return open_func(path, mode=mode, buffering=buffering, 
#                          encoding=encoding, errors=errors,
#                          newline=newline, closefd=closefd, 
#                          opener=opener)
#     return open_patched


# @pytest.fixture(autouse=True)
# def cleanup_files(monkeypatch):
#     files = []
#     monkeypatch.setattr(builtins, 'open', patch_open(builtins.open, files))
#     monkeypatch.setattr(io, 'open', patch_open(io.open, files))
#     yield
#     for file in files:
#         os.remove(file)

@pytest.fixture
def slurm_batch_unitjob(tmp_path):
    sim_input = {}
    calc_input = {
        'job_params': {},
        'slurm_params': True,
        'interactive': False,
    }
    job = UnitJob(sim_input, calc_input, workdir=tmp_path)
    return job

def test_update_and_read_output(slurm_batch_unitjob):
    output_update = {'test_val': 1}
    slurm_batch_unitjob.mkworkdir()
    slurm_batch_unitjob.update_output(output_update)
    print(slurm_batch_unitjob.get_output())
    assert slurm_batch_unitjob.get_output().get('test_val', False) == 1
