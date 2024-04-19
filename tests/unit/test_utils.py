'''
Tests for utils.py
'''
from pathlib import Path
import os
import pytest
from ase.io import read
from asimtools.utils import (
    join_names,
    get_atoms,
    get_images,
    change_dict_value,
    change_dict_values,
    parse_slice,
    write_yaml,
    read_yaml,
    get_env_input,
    get_calc_input,
    find_nth,
    get_nth_label,
    get_str_btn,
)
import ase.build

# Path to test data
STRUCT_DIR = Path(os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    '../data/structures/',
))

@pytest.mark.parametrize("test_input, expected",[
    (['a', 'b'],'a__b__'),
    (['a', 'b', 'c'],'a__b__c__'),
    (['_a', 'b.c'],'a__b.c__'),
    (['a', '-b.--'],'a__b__'),
    ([' ', '-b', 'c '],'b__c__'),
    (['', 'b'],'b__'),
])
def test_join_names(test_input, expected):
    ''' Test join_names '''
    assert join_names(test_input) == expected

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
    ({'image_file': str(STRUCT_DIR / 'images.xyz'), 'index': ':2'},
     [ase.build.bulk('Ar'), ase.build.bulk('Cu')]),
    ({'image_file': str(STRUCT_DIR / 'Li1.xyz')},
     [ase.build.bulk('Li').repeat((1,1,1))]),
    ({'pattern': str(STRUCT_DIR / 'Li1*.xyz')},
     [read(str(STRUCT_DIR / 'Li1.xyz'))]),
    ({'pattern': str(STRUCT_DIR / 'Li*.xyz')},
     [read(str(STRUCT_DIR / 'Li1.xyz')),
      read(str(STRUCT_DIR / 'Li2.xyz')),
      read(str(STRUCT_DIR / 'Li3.xyz'))]),
    ({'pattern': str(STRUCT_DIR / '*Li*.xyz'), 'skip_failed': True},
     [read(str(STRUCT_DIR / 'Li1.xyz')),
      read(str(STRUCT_DIR / 'Li2.xyz')),
      read(str(STRUCT_DIR / 'Li3.xyz'))]),
    ({'patterns': [str(STRUCT_DIR / 'Li1.xyz'), str(STRUCT_DIR / 'images.xyz')]},
     [read(str(STRUCT_DIR / 'Li1.xyz')),
      ase.build.bulk('Ar'),
      ase.build.bulk('Cu'),
      ase.build.bulk('Fe')]),
    ({'patterns': [
        str(STRUCT_DIR / 'bad_Li.xyz'), str(STRUCT_DIR / 'images.xyz')
    ], 'skip_failed': True},
     [ase.build.bulk('Ar'), ase.build.bulk('Cu'), ase.build.bulk('Fe')]),
    ({'patterns': [
        str(STRUCT_DIR / 'Li1.xyz'), str(STRUCT_DIR / 'images.xyz')
    ], 'index': -1},
     [read(str(STRUCT_DIR / 'Li1.xyz')), ase.build.bulk('Fe')]),
    ({'images': [ase.build.bulk('Ar'), ase.build.bulk('Fe')]},
     [ase.build.bulk('Ar'), ase.build.bulk('Fe')]),
])
def test_get_images(test_input, expected):
    ''' Test getting iterable of atoms from different inputs '''
    print('++input:', get_images(**test_input))
    print('++expected:', expected)
    input_images = get_images(**test_input)
    assert len(input_images) == len(expected)
    for image in input_images:
        assert image in expected

@pytest.mark.parametrize("test_input, expected",[
    (['l1', 'l2', 'l3'], {'l1': {'l2': {'l3': 'new_value'}}}),
    (['l1', 'l2'], {'l1': {'l2': 'new_value'}}),
    (['l1'], {'l1': 'new_value'}),
])
def test_change_dict_value(test_input, expected):
    ''' Test getting iterable of atoms from different inputs '''
    d = {'l1': {'l2': {'l3': 'l3_value'}}}
    new_d = change_dict_value(
        d, 'new_value', test_input, return_copy=True
    )
    assert new_d == expected
    assert new_d != d

@pytest.mark.parametrize("test_input, expected",[
    ([['l1', 'l21', 'l31'], ['l1', 'l21', 'l32']],
    {'l1': {'l21': {'l31': 'v1', 'l32': 'v2'}, 'l22': 'l22v'}}),
    ([['l1', 'l21', 'l31'], ['l1', 'l22']],
    {'l1': {'l21': {'l31': 'v1'}, 'l22': 'v2'}}),
    ([['l1', 'l21'], ['l1', 'l22']],
    {'l1': {'l21': 'v1', 'l22': 'v2'}}),
])
def test_change_dict_values(test_input, expected):
    ''' Test getting iterable of atoms from different inputs '''
    d = {'l1': {'l21': {'l31': 'l31v'}, 'l22': 'l22v'}}
    new_d = change_dict_values(
        d, ['v1', 'v2'], test_input, return_copy=True
    )
    assert new_d == expected
    assert new_d != d

@pytest.mark.parametrize("test_input, expected",[
    (':', slice(None, None, None)),
    ('1:', slice(1, None, None)),
    (':4', slice(None, 4, None)),
    ('1:5', slice(1, 5, None)),
    ('::-1', slice(None, None, -1)),
    ('3:67:4', slice(3, 67, 4)),
])
def test_parse_slice(test_input, expected):
    ''' Test parsing a slice '''
    assert parse_slice(test_input) == expected

def test_write_and_read_yaml(tmp_path):
    d = {'test_key': 'test_value'}
    yaml_file = tmp_path / 'test.yaml'
    write_yaml(yaml_file, d)
    assert yaml_file.exists()
    read_d = read_yaml(yaml_file)
    assert read_d == d

def test_get_env_input(tmp_path):
    cur_env_input = os.environ.get('ASIMTOOLS_ENV_INPUT', False)
    test_env_input = tmp_path / 'test_env_input.yaml'
    d = {'test_env': {'mode': {'slurm': False, 'interactive': True}}}
    write_yaml(test_env_input, d)
    os.environ['ASIMTOOLS_ENV_INPUT'] = str(test_env_input)
    env_input = get_env_input()
    if cur_env_input:
        os.environ['ASIMTOOLS_ENV_INPUT'] = cur_env_input
    assert env_input == d

def test_get_calc_input(tmp_path):
    cur_calc_input = os.environ.get('ASIMTOOLS_CALC_INPUT', False)
    test_calc_input = tmp_path / 'test_calc_input.yaml'
    d = {'test_calc': {
        'name': 'EMT', 'module': 'ase.calculators.emt', 'args': {}
    }}
    write_yaml(test_calc_input, d)
    os.environ['ASIMTOOLS_CALC_INPUT'] = str(test_calc_input)
    calc_input = get_calc_input()
    if cur_calc_input:
        os.environ['ASIMTOOLS_CALC_INPUT'] = cur_calc_input
    assert calc_input == d

@pytest.mark.parametrize("test_input, expected",[
    (('s','y'), 'ill'),
    (('illy','ally'), 's'),
    (('illy', 'l', 1), 'da'),
    (('ll', 'da', 1, 6), 'y'),
    ((None,'ll'), 'si'),
    (('ly',None,2), 'dally'),
])
def test_get_str_btn(test_input, expected):
    ''' Get string between two substrings '''
    s = 'sillysallywoulddillydally'
    assert get_str_btn(s, *test_input) == expected

@pytest.mark.parametrize("test_input, expected",[
    (1, 0),
    (2, 7),
    (3, 13),
])
def test_find_nth(test_input, expected):
    ''' Find nth occurence of substring in string '''
    s = 'measdfamesdfamebtnrtyh'
    assert find_nth(s, 'me', test_input) == expected

@pytest.mark.parametrize("test_input, expected",[
    (2, 'label-3'),
    (1, 'label2'),
    (0, 'label1'),
])
def test_get_nth_label(test_input, expected):
    ''' Find nth label in string '''
    s = '/path/to/array/id-00__label1__/id-0045__label2__/id-0023__label-3__/'
    assert get_nth_label(s, test_input) == expected

