'''
Tests for utils.py
'''
from pathlib import Path
import os
import pytest
from asimtools.utils import (
    join_names,
    write_yaml,
    get_atoms,
    get_images,
    parse_command_line
)
import ase.build

# Path to test data
STRUCT_DIR = Path(os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    '../data/structures/',
))

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

@pytest.mark.parametrize("test_input, expected",[
    ({'name': 'Ar'},
     ase.build.bulk('Ar')),
    ({'name': 'Ar', 'repeat': [2,2,2]},
     ase.build.bulk('Ar').repeat((2,2,2))),
    ({'name': 'H2O', 'builder': 'molecule'},
     ase.build.molecule('H2O')),
    ({'name': 'Li', 'crystalstructure': 'bcc', 'cubic': True},
     ase.build.bulk('Li', 'bcc', cubic=True)),
    ({'name': 'Kr', 'crystalstructure': 'bcc', 'a': 3},
     ase.build.bulk('Kr', 'bcc', a=3)),
    ({'symbol': 'Si', 'n': 6, 'm': 0, 'length': 4, 'builder': 'nanotube'},
     ase.build.nanotube(6, 0, length=4, symbol='Si')),
    ({'symbol': 'Al', 'size': (2,2,3), 'builder': 'fcc111'},
     ase.build.fcc111('Al', size=(2,2,3))),
    ({'image_file': STRUCT_DIR / 'Ar.xyz'},
     ase.build.bulk('Ar')),
    ({'image_file': STRUCT_DIR / 'Ar', 'format': 'cfg'},
     ase.build.bulk('Ar')),
    ({'atoms': ase.build.bulk('Ar')}, ase.build.bulk('Ar')),
])
def test_get_atoms(test_input, expected):
    ''' Test getting atoms from different inputs '''
    assert get_atoms(**test_input) == expected

@pytest.mark.parametrize("test_input, expected",[
    ({'image_file': str(STRUCT_DIR / 'images.xyz')},
     [ase.build.bulk('Ar'), ase.build.bulk('Cu'), ase.build.bulk('Fe')]),
    ({'pattern': str(STRUCT_DIR / 'Li*.xyz')},
     [ase.build.bulk('Li').repeat((1,1,1)),
      ase.build.bulk('Li').repeat((2,2,2)),
      ase.build.bulk('Li').repeat((3,3,3))]),
    ({'images': [ase.build.bulk('Ar'), ase.build.bulk('Fe')]},
     [ase.build.bulk('Ar'), ase.build.bulk('Fe')]),
])
def test_get_images(test_input, expected):
    ''' Test getting iterable of atoms from different inputs '''
    print(get_images(**test_input), expected, '*********')
    assert get_images(**test_input) == expected
