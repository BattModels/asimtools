'''
Tests for utils.py
'''
from pathlib import Path
from unittest.mock import patch, MagicMock
import os
import pytest
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from ase.io import read
from ase.calculators.emt import EMT
from pymatgen.core import Structure, Molecule, IStructure, IMolecule
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
    expand_wildcards,
    write_atoms,
    repeat_to_n,
    strip_symbols,
    get_axis_lims,
    improve_plot,
    write_csv_from_dict,
    new_db,
    check_if_slurm_job_is_running,
    find_files_by_regex,
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
    ({'image_file': STRUCT_DIR / 'Ar.xyz'},
     ase.build.bulk('Ar')),
    ({'atoms': ase.build.bulk('Ar')}, ase.build.bulk('Ar')),
    ({  
        'interface': 'pymatgen',
        'builder': 'Structure',
        'lattice': [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
        'coords': [[0, 0, 0]],
        'species': ['Ar'],
        'return_type': 'pymatgen',
    }, Structure(
        lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 3]],
        coords=[[0, 0, 0]],
        species=['Ar'],
        coords_are_cartesian=False,
    )),
    ({  
        'interface': 'pymatgen',
        'builder': 'Structure',
        'lattice': {'a': 3, 'b': 4, 'c': 5, 'alpha': 90, 'beta': 90, 'gamma': 90},
        'coords': [[0, 0, 0]],
        'species': ['Ar'],
        'return_type': 'pymatgen',
    }, Structure(
        lattice=[[3, 0, 0], [0, 4, 0], [0, 0, 5]],
        coords=[[0, 0, 0]],
        species=['Ar'],
        coords_are_cartesian=False,
    )),
    ({  
        'interface': 'pymatgen',
        'builder': 'IStructure',
        'lattice': {'a': 3, 'b': 4, 'c': 5, 'alpha': 90, 'beta': 90, 'gamma': 90},
        'coords': [[0, 0, 0]],
        'species': ['Ar'],
        'return_type': 'pymatgen',
    }, IStructure(
        lattice=[[3, 0, 0], [0, 4, 0], [0, 0, 5]],
        coords=[[0, 0, 0]],
        species=['Ar'],
        coords_are_cartesian=False,
    )),
    ({
        'interface': 'pymatgen',
        'builder': 'Molecule',
        'coords': [[0, 0, 0], [1.5, 1.5, 1.5], [1.5, -1.5, -1.5]],
        'species': ['O', 'H', 'H'],
        'spin_multiplicity': 1,
        'return_type': 'pymatgen',
    }, Molecule(
        species=['O', 'H', 'H'],
        coords=[[0, 0, 0], [1.5, 1.5, 1.5], [1.5, -1.5, -1.5]],
        spin_multiplicity=1,
    )),
    ({
        'interface': 'pymatgen',
        'builder': 'IMolecule',
        'coords': [[0, 0, 0], [1.5, 1.5, 1.5], [1.5, -1.5, -1.5]],
        'species': ['O', 'H', 'H'],
        'spin_multiplicity': 1,
        'return_type': 'pymatgen',
    }, IMolecule(
        species=['O', 'H', 'H'],
        coords=[[0, 0, 0], [1.5, 1.5, 1.5], [1.5, -1.5, -1.5]],
        spin_multiplicity=1,
    )),

])
def test_get_atoms(test_input, expected):
    ''' Test getting atoms from different inputs '''
    print('e', expected)
    assert get_atoms(**test_input) == expected

def test_get_atoms_constraints(tmp_path):
    atoms = ase.build.bulk('Cu').repeat((2,2,2))
    constrained_atoms = get_atoms(
        atoms=atoms, 
        constraints=[{'constraint': 'FixAtoms', 'indices': [0, 1]}]
    )
    assert len(constrained_atoms.constraints) == 1

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
    input_images = get_images(**test_input)
    assert len(input_images) == len(expected)
    for image in input_images:
        assert image in expected

@pytest.mark.parametrize("test_input, expected",[
    (
        {'image_file': str(STRUCT_DIR / 'adslab.xyz')},
        ('surface_properties', 'bulk_wyckoff'),
    ),
])
def test_write_atoms(test_input, expected, tmp_path):
    init_atoms = get_atoms(**test_input)
    testfile = tmp_path / 'test.xyz'
    write_atoms(testfile, init_atoms)

    with open(testfile, 'r') as f:
        lines = f.readlines()
    
    for prop in expected:
        assert prop in lines[1], f'"{prop}" not in file header'

def test_write_atoms_constraints(tmp_path):
    atoms = ase.build.bulk('Cu').repeat((2,2,2))
    atoms.set_constraint(ase.constraints.FixAtoms(indices=[0, 1]))
    write_atoms(tmp_path / 'test_constraints.xyz', atoms)
    read_atoms = read(tmp_path / 'test_constraints.xyz')
    assert len(read_atoms.constraints) > 0

@pytest.mark.parametrize("test_input, expected",[
    (
        {'name': 'Cu', 'interface': 'ase', 'builder': 'bulk'},
        ('forces', 'energy', 'stress'),
    ),
])
def test_write_atoms_calc(test_input, expected, tmp_path):
    init_atoms = get_atoms(**test_input)
    init_atoms.calc = EMT()
    init_atoms.get_potential_energy()
    testfile = tmp_path / 'test.xyz'
    write_atoms(testfile, init_atoms)

    with open(testfile, 'r') as f:
        lines = f.readlines()
    
    for prop in expected:
        assert prop in lines[1], f'"{prop}" not in file header'

def test_write_atoms_magmoms(tmp_path):
    ''' Test write_atoms. 
    Also test that when magmoms are provided by ASE which makes them nans
    and kills jobs using VASP etc., we change them to zero '''
    slab_ = read(STRUCT_DIR / 'adslab.xyz')
    write_atoms(tmp_path / 'slab.xyz', slab_)
    slab = read(tmp_path / 'slab.xyz')
    assert not np.any(np.isnan(slab.arrays['initial_magmoms']))
    assert 'surface_properties' in slab.arrays

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
    (['l1', 'l2', 'l3'], {'l1': {'l2': {'l3': 'l3_NEW'}}}),
])
def test_change_dict_value_placeholder(test_input, expected):
    ''' Test getting iterable of atoms from different inputs '''
    d = {'l1': {'l2': {'l3': 'l3_PLACEHOLDER'}}}
    new_d = change_dict_value(
        d, 'NEW', test_input, return_copy=True, placeholder='PLACEHOLDER',
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
    ((':', ), slice(None, None, None)),
    (('1:', ), slice(1, None, None)),
    ((':4', ), slice(None, 4, None)),
    (('1:5', ), slice(1, 5, None)),
    (('::-1', ), slice(None, None, -1)),
    (('3:67:4', ), slice(3, 67, 4)),
    (('1:', True), '$(seq 1 1 $END)'),
    (('1:5', True), '$(seq 1 1 5)'),
    (('::-1', True), '$(seq 0 -1 $END)'),
    (('3:67:4', True), '$(seq 3 4 67)'),
])
def test_parse_slice(test_input, expected):
    ''' Test parsing a slice '''
    assert parse_slice(*test_input) == expected

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

@pytest.mark.parametrize("test_input, expected",[
    ({'a': 'asdf'}, {'a': 'asdf'}),
    ({'a': 'asd*'}, {'a': 'asdf'}),
    ({'a': {'b': 'asd*'}}, {'a': {'b': 'asdf'}}),
    ({'a': {'b': 'asd*', 'c': 'zx*'}}, {'a': {'b': 'asdf', 'c': 'zxcv'}}),
])
def test_expand_wildcards(test_input, expected, tmp_path):
    ''' Test expanding paths '''
    fls = ['asdf', 'ascf', 'zxcv']
    for fl in fls:
        with open(tmp_path / fl, 'w') as f:
            f.write('')
    print(f'Found paths in {os.getcwd()}: {[f for f in Path(tmp_path).glob("*")]}')
    
    assert expand_wildcards(test_input, root_path=tmp_path) == expected

def test_repeat_to_n():
    ''' Test repeating unit cell to at least N atoms '''
    atoms = ase.build.bulk('Cu', crystalstructure='fcc', cubic=True, a=2.0)
    repeated_atoms = repeat_to_n(atoms, 16)
    assert len(repeated_atoms) == 16
    assert np.abs(repeated_atoms.get_cell()[0][0] - 2*2.0) < 1e-6
    assert np.abs(repeated_atoms.get_cell()[1][1] - 2*2.0) < 1e-6
    assert np.abs(repeated_atoms.get_cell()[2][2] - 1*2.0) < 1e-6
    assert len(repeat_to_n(atoms, 15)) == 16
    with pytest.raises(ValueError):
        repeat_to_n(atoms, 16, max_dim=4)


@pytest.mark.parametrize("test_input, expected", [
    ('hello', 'hello'),
    ('_hello_', 'hello'),
    ('-hello-', 'hello'),
    ('.hello.', 'hello'),
    (' hello ', 'hello'),
    ('_-hello-_', 'hello'),
    ('', ''),
])
def test_strip_symbols(test_input, expected):
    ''' Test stripping leading/trailing bad symbols from a string '''
    assert strip_symbols(test_input) == expected


@pytest.mark.parametrize("x, y, padding, expected_min, expected_max", [
    ([0, 1], [0, 1], 0.1, -0.1, 1.1),
    ([0, 2], [1, 3], 0.0, 0.0, 3.0),
    ([-1, 1], [-1, 1], 0.5, -2.0, 2.0),
])
def test_get_axis_lims(x, y, padding, expected_min, expected_max):
    ''' Test axis limit computation with padding '''
    lims = get_axis_lims(x, y, padding=padding)
    assert np.isclose(lims[0], expected_min)
    assert np.isclose(lims[1], expected_max)


def test_improve_plot():
    ''' Test that improve_plot sets label and tick fontsizes '''
    fig, ax = plt.subplots()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Title')
    fontsize = 16
    improve_plot(ax=ax, fontsize=fontsize)
    assert ax.xaxis.label.get_fontsize() == fontsize + 2
    assert ax.yaxis.label.get_fontsize() == fontsize + 2
    assert ax.title.get_fontsize() == fontsize + 4
    plt.close(fig)


def test_improve_plot_with_legend():
    ''' Test that improve_plot handles legends '''
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], label='line')
    ax.legend(title='legend_title')
    improve_plot(ax=ax, fontsize=12)
    leg = ax.get_legend()
    assert leg is not None
    plt.close(fig)


@pytest.mark.parametrize("data, columns, header", [
    ({'a': [1, 2, 3], 'b': [4, 5, 6]}, None, ''),
    ({'a': [1, 2], 'b': [3, 4], 'c': [5, 6]}, ['a', 'c'], 'my header'),
])
def test_write_csv_from_dict(data, columns, header, tmp_path):
    ''' Test writing a dict to CSV and reading it back '''
    fpath = tmp_path / 'test.csv'
    result_df = write_csv_from_dict(fpath, data, columns=columns, header=header)

    used_columns = columns if columns is not None else list(data.keys())
    assert list(result_df.columns) == used_columns

    read_df = pd.read_csv(fpath, comment='#')
    assert list(read_df.columns) == used_columns
    for col in used_columns:
        assert list(read_df[col]) == data[col]

    if header:
        with open(fpath, 'r', encoding='utf-8') as f:
            first_line = f.readline()
        assert header in first_line


def test_new_db(tmp_path):
    ''' Test creating a new ASE database '''
    import ase.db
    dbpath = tmp_path / 'test.db'
    db = new_db(str(dbpath))
    assert dbpath.exists()
    atoms = ase.build.bulk('Cu')
    db.write(atoms)
    assert db.count() == 1


def test_check_if_slurm_job_is_running_true():
    ''' Test that a job ID found in squeue stdout returns True '''
    mock_result = MagicMock()
    mock_result.stdout = ' 12345 some other text'
    with patch('asimtools.utils.subprocess.run', return_value=mock_result):
        assert check_if_slurm_job_is_running(12345) is True


def test_check_if_slurm_job_is_running_false():
    ''' Test that a job ID not in squeue stdout returns False '''
    mock_result = MagicMock()
    mock_result.stdout = ' 99999 some other text'
    with patch('asimtools.utils.subprocess.run', return_value=mock_result):
        assert check_if_slurm_job_is_running(12345) is False


# ── find_files_by_regex ───────────────────────────────────────────────────────

def test_find_files_by_regex_matches(tmp_path):
    ''' Files whose full path matches the regex are returned '''
    (tmp_path / 'bulk_Na.xyz').write_text('')
    (tmp_path / 'bulk_Li.xyz').write_text('')
    (tmp_path / 'surface_Na.xyz').write_text('')
    matches = find_files_by_regex(r'bulk_.*\.xyz', search_dir=tmp_path)
    assert len(matches) == 2
    assert all('bulk_' in m for m in matches)


def test_find_files_by_regex_no_match(tmp_path):
    ''' Returns empty list when nothing matches '''
    (tmp_path / 'bulk_Na.xyz').write_text('')
    matches = find_files_by_regex(r'surface_.*\.xyz', search_dir=tmp_path)
    assert matches == []


def test_find_files_by_regex_non_recursive(tmp_path):
    ''' Non-recursive search does not descend into subdirectories '''
    subdir = tmp_path / 'sub'
    subdir.mkdir()
    (subdir / 'bulk_Na.xyz').write_text('')
    matches = find_files_by_regex(r'bulk_.*\.xyz', search_dir=tmp_path, recursive=False)
    assert matches == []


def test_find_files_by_regex_recursive(tmp_path):
    ''' Recursive search finds files in subdirectories '''
    subdir = tmp_path / 'sub'
    subdir.mkdir()
    (subdir / 'bulk_Na.xyz').write_text('')
    matches = find_files_by_regex(r'bulk_.*\.xyz', search_dir=tmp_path, recursive=True)
    assert len(matches) == 1


def test_find_files_by_regex_natsorted(tmp_path):
    ''' Results are returned in natural sort order '''
    for name in ['bulk_10.xyz', 'bulk_2.xyz', 'bulk_1.xyz']:
        (tmp_path / name).write_text('')
    matches = find_files_by_regex(r'bulk_\d+\.xyz', search_dir=tmp_path)
    basenames = [Path(m).name for m in matches]
    assert basenames == ['bulk_1.xyz', 'bulk_2.xyz', 'bulk_10.xyz']


# ── get_atoms with regex ──────────────────────────────────────────────────────

def test_get_atoms_regex(tmp_path):
    ''' get_atoms regex finds a single matching file and reads it '''
    expected = ase.build.bulk('Ar')
    expected.write(str(tmp_path / 'bulk_Ar.xyz'), format='extxyz')
    atoms = get_atoms(regex=r'bulk_Ar\.xyz', regex_kwargs={'search_dir': tmp_path})
    assert atoms == expected


def test_get_atoms_regex_multiple_matches_raises(tmp_path):
    ''' get_atoms raises ValueError when regex matches more than one file '''
    for name in ['bulk_Ar.xyz', 'bulk_Cu.xyz']:
        ase.build.bulk('Ar').write(str(tmp_path / name), format='extxyz')
    with pytest.raises(ValueError, match='matched 2 files'):
        get_atoms(regex=r'bulk_.*\.xyz', regex_kwargs={'search_dir': tmp_path})


def test_get_atoms_regex_no_match_raises(tmp_path):
    ''' get_atoms raises ValueError when regex matches no files '''
    with pytest.raises(ValueError, match='matched 0 files'):
        get_atoms(regex=r'surface_.*\.xyz', regex_kwargs={'search_dir': tmp_path})


# ── get_images with regex ─────────────────────────────────────────────────────

def test_get_images_regex(tmp_path):
    ''' get_images regex finds all matching files and returns their atoms '''
    expected = [ase.build.bulk('Ar'), ase.build.bulk('Cu'), ase.build.bulk('Fe')]
    for atoms in expected:
        atoms.write(
            str(tmp_path / f'bulk_{atoms.get_chemical_formula()}.xyz'),
            format='extxyz',
        )
    result = get_images(regex=r'bulk_.*\.xyz', regex_kwargs={'search_dir': tmp_path})
    assert len(result) == 3
    for atoms in result:
        assert atoms in expected


def test_get_images_regex_skip_failed(tmp_path):
    ''' get_images regex skips unreadable files when skip_failed=True '''
    ase.build.bulk('Ar').write(str(tmp_path / 'bulk_Ar.xyz'), format='extxyz')
    (tmp_path / 'bulk_bad.xyz').write_text('not valid xyz content')
    result = get_images(
        regex=r'bulk_.*\.xyz',
        regex_kwargs={'search_dir': tmp_path},
        skip_failed=True,
    )
    assert len(result) == 1