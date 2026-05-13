from typing import Dict, Sequence, Optional, Union
from glob import glob
from pathlib import Path
from natsort import natsorted
import re
import numpy as np
from asimtools.utils import get_str_btn, get_nth_label, find_files_by_regex

def prepare_array_vals(
    key_sequence: Optional[Sequence[str]] = None,
    array_values: Optional[Sequence] = None,
    file_pattern: Optional[str] = None,
    file_regex: str | None = None,
    file_regex_kwargs: dict | None = None,
    linspace_args: Optional[Sequence] = None,
    arange_args: Optional[Sequence] = None,
    env_ids: Optional[Union[Sequence[str],str]] = None,
    labels: Optional[Union[Sequence,str]] = 'values',
    label_prefix: Optional[str] = None,
    str_btn_args: Optional[Dict] = None,
    regex_label_args: dict | None = None,
    secondary_key_sequences: Optional[Sequence] = None,
    secondary_array_values: Optional[Sequence] = None,
    as_integers: Optional[bool] = False,
):
    """Helper function for preparing things needed for the different arrays

    :param key_sequence: Sequence of keys to access value to be iterated over,
        defaults to None
    :type key_sequence: Optional[Sequence[str]], optional
    :param array_values: values to be iterated over in each simulation,
        defaults to None
    :type array_values: Optional[Sequence], optional
    :param file_pattern: glob pattern of files to be iterated over in each
        simulation, defaults to None
    :type file_pattern: Optional[str], optional
    :param file_regex: regex pattern matched against full file paths to be
        iterated over in each simulation, defaults to None
    :type file_regex: str, optional
    :param file_regex_kwargs: kwargs forwarded to
        :func:`asimtools.utils.find_files_by_regex` (e.g. ``search_dir``,
        ``recursive``), defaults to None
    :type file_regex_kwargs: dict, optional
    :param linspace_args: arguments to pass to :func:`numpy.linspace` to be
        iterated over in each simulation, defaults to None
    :type linspace_args: Optional[Sequence], optional
    :param arange_args: arguments to pass to :func:`numpy.arange` to be
        iterated over in each simulation, defaults to None
    :type arange_args: Optional[Sequence], optional
    :param labels: Custom labels to use for each simulation. If ``'str_btn'``,
        provide arguments via ``str_btn_args``. If ``'regex'``, provide a
        pattern and optional group via ``regex_label_args``. If an integer N,
        the Nth label in the file_pattern or array_values is used, defaults to
        ``'values'``
    :type labels: Sequence, optional
    :param regex_label_args: Dict with keys ``pattern`` (required) and
        ``group`` (optional, default 1) used when ``labels='regex'`` to extract
        a label from each array value via :func:`re.search`, defaults to None
    :type regex_label_args: dict, optional
    :param label_prefix: Prefix to add before labels which can make extracting 
        data from file paths easier, defaults to None
    :type label_prefix: str, optional
    :param env_ids: Environment(s) to be used for each simulation, must either
        be a list with as many env_ids as array values or a string with the
        env_id to be used by all simulations, defaults to None
    :type env_ids: Optional[Union[Sequence[str],str]], optional
    :param secondary_key_sequences: list of other keys to iterate over in
        tandem with key_sequence to allow changing multiple key-value pairs,
        defaults to None
    :type secondary_key_sequences: Sequence, optional
    :param secondary_array_values: list of other other array_values to iterate
        over in tandem with array_values to allow changing multiple key-value
        pairs, defaults to None
    :type secondary_array_values: Sequence, optional
    :param as_integers: Whether to return the values as integers, useful for 
        indexing
    :type as_integers: bool, False
    :return: Results
    :rtype: Dict
    """
    possible_values = [
        array_values, file_pattern, file_regex, linspace_args, arange_args
    ]
    possible_strs = [
        'array_values', 'file_pattern', 'file_regex', 'linspace_args', 'arange_args'
    ]
    assert possible_values.count(None) + 1 == len(possible_values), \
        f'Provide only one of {[str(v) for v in possible_strs]}'

    if file_pattern is not None:
        array_values = natsorted(glob(str(file_pattern)))
        assert len(array_values) > 0, f'No file_pattern matching {file_pattern}'
    elif file_regex is not None:
        array_values = find_files_by_regex(file_regex, **(file_regex_kwargs or {}))
        assert len(array_values) > 0, f'No files matching regex "{file_regex}"'
    elif linspace_args is not None:
        array_values = np.linspace(*linspace_args)
        array_values = [float(v) for v in array_values]
    elif arange_args is not None:
        array_values = np.arange(*arange_args)
        if as_integers:
            array_values = [int(v) for v in array_values]
        else:
            array_values = [float(v) for v in array_values]

    assert len(array_values) > 0, f'No array_values found'

    if labels == 'str_btn':
        assert str_btn_args is not None, 'Provide str_btn_args for labels'
        labels = [get_str_btn(s, *str_btn_args) for s in array_values]
    elif labels == 'regex':
        assert regex_label_args is not None, 'Provide regex_label_args for labels'
        pattern = regex_label_args['pattern']
        group = regex_label_args.get('group', 1)
        labels = []
        for v in array_values:
            m = re.search(pattern, str(v))
            if m is None:
                raise ValueError(
                    f'regex pattern "{pattern}" did not match value "{v}"'
                )
            labels.append(m.group(group))
    elif labels == 'values':
        if key_sequence is not None:
            labels = [f'{key_sequence[-1]}-{val}' for val in array_values]
        else:
            labels = [f'value-{val}' for val in array_values]
    elif isinstance(labels, int):
        labels = [get_nth_label(s, labels) for s in array_values]
    elif labels is None:
        labels = [str(i) for i in range(len(array_values))]

    # In case the label is a file path with /  or space characters
    labels = [label.replace('/', '+') for label in labels]
    labels = [label.replace(' ', '+') for label in labels]
    labels = [label.replace('*', '') for label in labels]

    if label_prefix is not None:
        labels = [label_prefix + '-' + label for label in labels]

    assert len(labels) == len(array_values), \
        f'Num. of array_values ({len(array_values)}) must match num.'\
        f'of labels ({len(labels)})'
    
    # File patterns should be resolved fully for placeholders to work
    if file_pattern is not None or file_regex is not None:
        array_values = [str(Path(v).resolve()) for v in array_values]

    if secondary_array_values is not None:
        nvals = len(secondary_array_values)
        nkeys = len(secondary_key_sequences)
        assert nvals == nkeys, \
            f'Num. of secondary values ({nvals}) does not match num. of '\
            f'secondary keys ({nkeys})'
        for l in secondary_array_values:
            assert len(l) == len(labels), \
                f"Secondary values ({len(l)}) not same length as array values"\
                f" ({len(labels)})"

    if env_ids is None:
        pass
    elif isinstance(env_ids, str):
        env_ids = [env_ids for i in range(len(labels))]
    else:
        assert len(env_ids) == len(labels) or isinstance(env_ids, str), \
            'Provide one env_id or as many as there are array_values'

    results = {
        'array_values': array_values,
        'labels': labels,
        'env_ids': env_ids,
        'secondary_array_values': secondary_array_values,
        'secondary_key_sequences': secondary_key_sequences,
    }

    return results
