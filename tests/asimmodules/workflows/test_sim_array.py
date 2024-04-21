"""
Tests for running asimmodules using asim_run.py
"""
import os
from pathlib import Path
import pytest
from asimtools.job import load_job_from_directory
from asimtools.asimmodules.workflows.sim_array import sim_array

# Path to test data
FILE_DIR = Path(os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    '../../data/arrays/',
))

@pytest.mark.parametrize("sim_input",["do_nothing_sim_input"])
@pytest.mark.parametrize("env_input",["all_env_input"])
@pytest.mark.parametrize(
    "array_src, expected_dirs", [
        (
            {
                "array_values": [0,1,2,3]
            },
            [
                'id-0000__value-0__',
                'id-0001__value-1__',
                'id-0002__value-2__',
                'id-0003__value-3__'
            ]
        ),
        (
            {
                "array_values": ['0','1','2','3'],
                "labels": ['a', 'b', 'c', 'd'],
            },
            ['id-0000__a__','id-0001__b__','id-0002__c__','id-0003__d__']
        ),
        (
            {
                "file_pattern": str(FILE_DIR / 'test_file_*'),
                "labels": "str_btn",
                "str_btn_args": ("test_", None)
            },
            [
                'id-0000__file_0__',
                'id-0001__file_1__',
                'id-0002__file_2__',
                'id-0003__file_3__',
            ]
        ),
        (
            {
                "arange_args": (0,4,1)
            },
            [
                'id-0000__value-0.0__',
                'id-0001__value-1.0__',
                'id-0002__value-2.0__',
                'id-0003__value-3.0__'
            ]
        ),
        (
            {
                "linspace_args": (0,3,4)
            },
            [
                'id-0000__value-0.0__',
                'id-0001__value-1.0__',
                'id-0002__value-2.0__',
                'id-0003__value-3.0__'
            ]
        ),
    ]
)
@pytest.mark.parametrize(
    "test_env_ids, expected_env_ids", [
        ("inline", ["inline"] * 4),
        (
            ["inline","inline2","inline","inline2"],
            ["inline","inline2","inline","inline2"]
        )
    ]
)
@pytest.mark.parametrize(
    "test_secondary_args", [
        {
            "secondary_key_sequences": [["args","duration"]],
            "secondary_array_values": [[0.1,0.2,0.3,0.4]],
        },
        {
            "secondary_key_sequences": [["job_name"]],
            "secondary_array_values": [["one","two","three","four"]],
        },
    ]
)
def test_sim_array(
    sim_input,
    env_input,
    array_src,
    expected_dirs,
    test_env_ids,
    expected_env_ids,
    test_secondary_args,
    tmp_path,
    request,
):
    cur_dir = Path('.').resolve()
    wdir = tmp_path
    os.chdir(str(wdir))
    sim_input = request.getfixturevalue(sim_input)
    sim_input["args"]["duration"] = 0.1
    env_input = request.getfixturevalue(env_input)
    sim_array(
        template_sim_input=sim_input,
        env_input=env_input,
        env_ids=test_env_ids,
        **array_src,
        **test_secondary_args,
    )

    for i in range(4):
        unitjob_workdir = Path(wdir / expected_dirs[i]).resolve()
        assert unitjob_workdir.exists()
        assert (unitjob_workdir / "sim_input.yaml").exists()
        assert (unitjob_workdir / "calc_input.yaml").exists()
        assert (unitjob_workdir / "env_input.yaml").exists()

        unitjob = load_job_from_directory(unitjob_workdir)
        uj_sim_input = unitjob.get_sim_input()
        uj_output = unitjob.get_output()

        assert uj_output["status"] == "complete"
        assert int(uj_sim_input["distributed_id"]) == i
        assert uj_sim_input['env_id'] == expected_env_ids[i]

        for j, _ in enumerate(test_secondary_args['secondary_key_sequences']):
            val = uj_sim_input.copy()
            for k in test_secondary_args['secondary_key_sequences'][j]:
                val = val[k]
            assert val == test_secondary_args['secondary_array_values'][j][i]

    os.chdir(str(cur_dir))
