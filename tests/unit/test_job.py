'''
Test Job class
'''
#pylint: disable=missing-function-docstring
#pylint: disable=redefined-outer-name

from asimtools.job import UnitJob

def test_gen_input_files(
    slurm_batch_calc_input,
    singlepoint_sim_input,
    tmp_path
):
    wdir = tmp_path / 'wdir'
    unitjob = UnitJob(
        slurm_batch_calc_input,
        singlepoint_sim_input,
        wdir,
    )
    unitjob.gen_input_files()
    assert wdir.exists()
    assert (wdir / 'sim_input.yaml').exists()
    assert (wdir / 'calc_input.yaml').exists()
    assert (wdir / 'test_output.yaml').exists()
    assert (wdir / 'job.sh').exists()


def test_update_and_read_output(
    interactive_calc_input,
    singlepoint_sim_input,
    tmp_path
):
    output_update = {'test_val': 1}
    unitjob = UnitJob(
        interactive_calc_input,
        singlepoint_sim_input,
        tmp_path,
    )
    unitjob.gen_input_files()
    unitjob.update_output(output_update)
    assert unitjob.get_output().get('test_val', False) == 1
