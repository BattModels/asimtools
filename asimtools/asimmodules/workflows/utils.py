from typing import Dict, Sequence, Optional, Union
from glob import glob
import numpy as np
from asimtools.utils import get_str_btn

def prepare_array_vals(
    key_sequence: Optional[Sequence[str]] = None,
    array_values: Optional[Sequence] = None,
    file_pattern: Optional[str] = None,
    linspace_args: Optional[Sequence] = None,
    arange_args: Optional[Sequence] = None,
    env_ids: Optional[Union[Sequence[str],str]] = None,
    labels: Optional[Union[Sequence,str]] = 'values',
    label_prefix: Optional[str] = None,
    str_btn_args: Optional[Dict] = None,
    secondary_key_sequences: Optional[Sequence] = None,
    secondary_array_values: Optional[Sequence] = None,
):
    """Helper function for preparing things needed for the different arrays

    :param key_sequence: Sequence of keys to access value to be iterated over,
        defaults to None
    :type key_sequence: Optional[Sequence[str]], optional
    :param array_values: values to be iterated over in each simulation,
        defaults to None
    :type array_values: Optional[Sequence], optional
    :param file_pattern: pattern of files to be iterated over in each
        simulation, defaults to None
    :type file_pattern: Optional[str], optional
    :param linspace_args: arguments to pass to :func:`numpy.linspace` to be
        iterated over in each simulation, defaults to None
    :type linspace_args: Optional[Sequence], optional
    :param arange_args: arguments to pass to :func:`numpy.arange` to be
        iterated over in each simulation, defaults to None
    :type arange_args: Optional[Sequence], optional
    :param labels: Custom labels to use for each simulation, defaults to None
    :type labels: Sequence, optional
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
    :return: Results
    :rtype: Dict
    """
    possible_values = [
        array_values, file_pattern, linspace_args, arange_args
    ]
    possible_strs = [
        'array_values','file_pattern','linspace_args','arange_args'
    ]
    assert possible_values.count(None) + 1 == len(possible_values), \
        f'Provide only one of {[str(v) for v in possible_strs]}'

    if file_pattern is not None:
        array_values = sorted(glob(str(file_pattern)))
    elif linspace_args is not None:
        array_values = np.linspace(*linspace_args)
        array_values = [float(v) for v in array_values]
    elif arange_args is not None:
        array_values = np.arange(*arange_args)
        array_values = [float(v) for v in array_values]

    assert len(array_values) > 0, 'No array values or files found'

    if labels == 'str_btn':
        assert str_btn_args is not None, 'Provide str_btn_args for labels'
        labels = [get_str_btn(s, *str_btn_args) for s in array_values]
    elif labels == 'values':
        if key_sequence is not None:
            labels = [f'{key_sequence[-1]}-{val}' for val in array_values]
        else:
            labels = [f'value-{val}' for val in array_values]
    elif labels is None:
        labels = [str(i) for i in range(len(array_values))]

    # In case the label is a file path with / characters
    labels = [label.replace('/', '+') for label in labels]

    if label_prefix is not None:
        labels = [label_prefix + '-' + label for label in labels]

    assert len(labels) == len(array_values), \
        'Num. of array_values must match num. of labels'

    if secondary_array_values is not None:
        nvals = len(secondary_array_values)
        nkeys = len(secondary_key_sequences)
        assert nvals == nkeys, \
            f'{nvals} secondary values does not match {nkeys} secondary keys'
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
