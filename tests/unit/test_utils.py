'''
Tests for utils.py
'''
import pytest
from asimtools.utils import join_names, write_yaml
from asimtools.scripts.asim_run import parse_command_line

@pytest.mark.parametrize("test_input, expected",[
    (['a', 'b'],'a_b'),
    (['a', 'b', 'c'],'a_b_c'),
    (['_a', 'b.c'],'a_b.c'),
    (['a', '-b.--'],'a_b'),
    ([' ', '-b', 'c '],'b_c'),
    (['', 'b'],'b'),
])
def test_join_names(test_input, expected):
    ''' Test join_names '''
    assert join_names(test_input) == expected

def test_parse_command_line(
    tmp_path,
    lj_interactive_calc_input,
    singlepoint_sim_input
):
    ''' Tests commandline argument parsing '''
    calc_input_file = str(tmp_path / 'calc_input.yaml')
    sim_input_file = str(tmp_path / 'sim_input.yaml')
    write_yaml(calc_input_file, lj_interactive_calc_input)
    write_yaml(sim_input_file, singlepoint_sim_input)

    calc_config_parse, sim_config_parse = parse_command_line(
        [calc_input_file, sim_input_file]
    )

    ## some checks
    assert lj_interactive_calc_input == calc_config_parse
    assert singlepoint_sim_input == sim_config_parse

#TODO
def test_get_atoms():
    ''' Test getting atoms from different inputs '''
    pass

#TODO
def test_get_images():
    ''' Test getting iterable of atoms from different inputs '''
    pass
