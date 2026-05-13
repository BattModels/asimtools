'''
Tests for asimmodules/workflows/utils.py
'''
import os
from pathlib import Path
import pytest
from asimtools.asimmodules.workflows.utils import prepare_array_vals

FILE_DIR = Path(os.path.dirname(os.path.realpath(__file__))) / '../data/arrays'

# Pre-resolved absolute paths for the four test fixture files
_FILE_PATHS = [
    str((FILE_DIR / f'test_file_{i}__').resolve()) for i in range(4)
]


# ── happy-path parametrize ────────────────────────────────────────────────────

@pytest.mark.parametrize("kwargs, expected_values, expected_labels", [
    # array_values — default 'values' labels use key_sequence tail
    (
        {'array_values': [0, 1, 2], 'key_sequence': ['param']},
        [0, 1, 2],
        ['param-0', 'param-1', 'param-2'],
    ),
    # array_values — explicit label list
    (
        {'array_values': [0, 1, 2], 'labels': ['a', 'b', 'c']},
        [0, 1, 2],
        ['a', 'b', 'c'],
    ),
    # array_values — labels=None produces index strings
    (
        {'array_values': [0, 1, 2], 'labels': None},
        [0, 1, 2],
        ['0', '1', '2'],
    ),
    # linspace_args
    (
        {'linspace_args': (0, 3, 4)},
        [0.0, 1.0, 2.0, 3.0],
        ['value-0.0', 'value-1.0', 'value-2.0', 'value-3.0'],
    ),
    # arange_args — floats by default
    (
        {'arange_args': (0, 4, 1)},
        [0.0, 1.0, 2.0, 3.0],
        ['value-0.0', 'value-1.0', 'value-2.0', 'value-3.0'],
    ),
    # arange_args — as_integers=True
    (
        {'arange_args': (0, 4, 1), 'as_integers': True},
        [0, 1, 2, 3],
        ['value-0', 'value-1', 'value-2', 'value-3'],
    ),
    # file_pattern — labels derived via str_btn
    (
        {
            'file_pattern': str(FILE_DIR / 'test_file_*'),
            'labels': 'str_btn',
            'str_btn_args': ('test_', None),
        },
        _FILE_PATHS,
        ['file_0__', 'file_1__', 'file_2__', 'file_3__'],
    ),
    # file_regex — same files via regex, same labels
    (
        {
            'file_regex': r'test_file_\d',
            'file_regex_kwargs': {'search_dir': FILE_DIR},
            'labels': 'str_btn',
            'str_btn_args': ('test_', None),
        },
        _FILE_PATHS,
        ['file_0__', 'file_1__', 'file_2__', 'file_3__'],
    ),
])
def test_prepare_array_vals(kwargs, expected_values, expected_labels):
    ''' prepare_array_vals returns correct array_values and labels '''
    result = prepare_array_vals(**kwargs)
    assert result['array_values'] == expected_values
    assert result['labels'] == expected_labels


# ── env_ids ───────────────────────────────────────────────────────────────────

def test_prepare_array_vals_env_ids_broadcast():
    ''' A single env_id string is broadcast to match array length '''
    result = prepare_array_vals(array_values=[0, 1, 2], env_ids='inline')
    assert result['env_ids'] == ['inline', 'inline', 'inline']


def test_prepare_array_vals_env_ids_list():
    ''' A list of env_ids passes through unchanged '''
    env_ids = ['a', 'b', 'c']
    result = prepare_array_vals(array_values=[0, 1, 2], env_ids=env_ids)
    assert result['env_ids'] == env_ids


# ── label_prefix ──────────────────────────────────────────────────────────────

def test_prepare_array_vals_label_prefix():
    ''' label_prefix is prepended to every label '''
    result = prepare_array_vals(
        array_values=[0, 1], labels=['x', 'y'], label_prefix='run'
    )
    assert result['labels'] == ['run-x', 'run-y']


# ── secondary values ──────────────────────────────────────────────────────────

def test_prepare_array_vals_secondary_values():
    ''' secondary_array_values and secondary_key_sequences pass through '''
    sec_vals = [[10, 20, 30]]
    sec_keys = [['a', 'b']]
    result = prepare_array_vals(
        array_values=[0, 1, 2],
        secondary_array_values=sec_vals,
        secondary_key_sequences=sec_keys,
    )
    assert result['secondary_array_values'] == sec_vals
    assert result['secondary_key_sequences'] == sec_keys


# ── file path resolution ──────────────────────────────────────────────────────

def test_prepare_array_vals_file_paths_resolved():
    ''' File-based array_values are resolved to absolute paths '''
    result = prepare_array_vals(file_pattern=str(FILE_DIR / 'test_file_*'))
    assert all(Path(v).is_absolute() for v in result['array_values'])


def test_prepare_array_vals_file_regex_paths_resolved():
    ''' file_regex array_values are also resolved to absolute paths '''
    result = prepare_array_vals(
        file_regex=r'test_file_\d',
        file_regex_kwargs={'search_dir': FILE_DIR},
    )
    assert all(Path(v).is_absolute() for v in result['array_values'])


# ── file_regex with custom search_dir ─────────────────────────────────────────

def test_prepare_array_vals_file_regex_search_dir(tmp_path):
    ''' file_regex_kwargs search_dir scopes the search correctly '''
    for name in ['run_0.xyz', 'run_1.xyz', 'other.txt']:
        (tmp_path / name).write_text('')
    result = prepare_array_vals(
        file_regex=r'run_\d\.xyz',
        file_regex_kwargs={'search_dir': tmp_path},
    )
    assert len(result['array_values']) == 2


# ── error cases ───────────────────────────────────────────────────────────────

def test_prepare_array_vals_mutual_exclusion_raises():
    ''' Providing two sources raises AssertionError '''
    with pytest.raises(AssertionError):
        prepare_array_vals(array_values=[0, 1], linspace_args=(0, 1, 2))


def test_prepare_array_vals_file_regex_no_match_raises(tmp_path):
    ''' file_regex with no matching files raises AssertionError '''
    with pytest.raises(AssertionError, match='No files matching regex'):
        prepare_array_vals(
            file_regex=r'nonexistent_.*\.xyz',
            file_regex_kwargs={'search_dir': tmp_path},
        )


def test_prepare_array_vals_file_pattern_no_match_raises():
    ''' file_pattern with no matching files raises AssertionError '''
    with pytest.raises(AssertionError, match='No file_pattern'):
        prepare_array_vals(file_pattern='/nonexistent/path/to/files/test_*')


# ── labels='regex' ────────────────────────────────────────────────────────────

def test_prepare_array_vals_regex_labels_numbered_group():
    ''' labels=regex extracts a numbered capture group from each value '''
    result = prepare_array_vals(
        array_values=['bulk_Na_relaxed', 'bulk_Li_relaxed'],
        labels='regex',
        regex_label_args={'pattern': r'bulk_(\w+)_relaxed', 'group': 1},
    )
    assert result['labels'] == ['Na', 'Li']


def test_prepare_array_vals_regex_labels_named_group():
    ''' labels=regex works with a named capture group '''
    result = prepare_array_vals(
        array_values=['run_001_final', 'run_042_final'],
        labels='regex',
        regex_label_args={'pattern': r'run_(?P<id>\d+)_final', 'group': 'id'},
    )
    assert result['labels'] == ['001', '042']


def test_prepare_array_vals_regex_labels_no_match_raises():
    ''' labels=regex raises ValueError when the pattern does not match a value '''
    with pytest.raises(ValueError, match='did not match value'):
        prepare_array_vals(
            array_values=['bulk_Na', 'surface_Li'],
            labels='regex',
            regex_label_args={'pattern': r'bulk_(\w+)_relaxed', 'group': 1},
        )


def test_prepare_array_vals_regex_labels_no_args_raises():
    ''' labels=regex raises AssertionError when regex_label_args is not provided '''
    with pytest.raises(AssertionError, match='Provide regex_label_args'):
        prepare_array_vals(array_values=[0, 1], labels='regex')
