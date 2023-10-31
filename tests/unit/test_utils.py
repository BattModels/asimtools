'''
Tests for utils.py
'''
from pathlib import Path
import os
import pytest
from asimtools.utils import (
    join_names,
    get_atoms,
    get_images,
    change_dict_value,
    parse_slice,
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
    ({'image_file': str(STRUCT_DIR / 'Li1.xyz')},
     [ase.build.bulk('Li').repeat((1,1,1))]),
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

