'''
Utilities and helper functions for reading and writing data using set
standards
'''
import getopt
import sys
from typing import TypeVar, Iterable, Tuple
from glob import glob
import yaml
import pandas as pd
from ase.io import read
from ase.build import bulk

Atoms = TypeVar('Atoms')

def read_yaml(yaml_path: str) -> dict:
    ''' Read a yaml file'''
    with open(yaml_path, 'r', encoding='utf-8') as f:
        output = yaml.safe_load(f)
    return output

def write_yaml(yaml_path: str, yaml_dict: dict) -> None:
    ''' Write a yaml file '''
    with open(yaml_path, 'w', encoding='utf-8') as f:
        yaml.dump(yaml_dict, f)

def write_csv_from_dict(fname: str, data: dict, columns: Iterable = None) -> pd.DataFrame:
    ''' Write a csv from entries in a dictionary '''

    if columns is None:
        columns = list(data.keys())

    csv_data = pd.DataFrame(data, columns=columns)
    csv_data.to_csv(fname, index=False)
    return csv_data

def strip_symbols(substr: str) -> str:
    ''' 
    Helper for join_names. Strips symbols used for special purposes
    '''
    bad_symbols = ['_', '-', '.', ' ']
    if len(substr) == 0:
        return ''
    if substr[0] in bad_symbols or substr[-1] in bad_symbols:
        for bad_symbol in bad_symbols:
            substr = substr.strip(bad_symbol)
        return strip_symbols(substr)
    return substr

def join_names(substrs: Iterable[str]) -> str:
    ''' 
    Join multiple parts of a file/dir name using a standard format where
    individual parts of a file name are separated by one underscore with no
    hyphens, periods or underscores next to an underscore
    '''
    name = ''
    new_substrs = []
    for substr in substrs:
        new_substrs.append(strip_symbols(substr))

    final_substrs = []
    for substr in new_substrs:
        if len(substr) > 0:
            final_substrs.append(substr)
    name = '_'.join(final_substrs)
    return name

def get_atoms(
    symbol: str = None,
    crystal_structure: str = None,
    input_file: str = None,
    atoms: Atoms = None,
    **kwargs
) -> Atoms:
    '''
    Helper for returning an Atoms object by ASE specification, input
    file or directly using an atoms object
    '''
    assert (crystal_structure is not None and symbol is not None) or \
            input_file is not None or \
            atoms is not None, \
            "Please specify atoms or input structure and symbol or file"

    if input_file is not None:
        assert isinstance(input_file, str), "Input file must be a path to\
            an ASE readable structure"
        atoms = read(input_file)
    elif atoms is None:
        assert kwargs.get('a', False), \
            "Please specify at least lattice parameter a"
        atoms = bulk(symbol, crystal_structure, **kwargs)
    else:
        assert atoms is not None, 'Please give an input structure'

    return atoms

def get_images(
    input_file: str = None,
    pattern: str = None,
    images: Iterable[Atoms] = None,
    format: str = None,
    index: slice = slice(0,-1),
) -> list:
    '''
    Helper for returning an Atoms object by ASE specification, input
    file or directly using an atoms object
    '''
    assert (input_file is not None) or \
        (pattern is not None) or \
        (images is not None), \
        "Please specify file, pattern or iterable"

    if input_file is not None:
        images = read(input_file, index=index, format=format)
    elif pattern is not None:
        image_files = glob(pattern)
        images = [
            read(image_file, format=format) for image_file in image_files
        ]
    elif images is not None:
        # images = eval(f'images[{index}]')
        images = images[index]
    else:
        images = []

    assert len(images) > 0, 'No images found'

    return images

def parse_command_line(argv) -> Tuple[dict, dict]:
    ''' Parse command line inputs for simulation script '''
    usage = 'python singlepoint.py calc_params.yaml sim_params.yaml [-h]'

    assert len(argv) >= 2, usage
    try:
        calc_inputfile = argv[0]
        calc_params = read_yaml(calc_inputfile)
    except IndexError:
        print('Failed to load calc_inputfile\n', usage)
        raise

    try:
        sim_inputfile = argv[1]
        sim_params = read_yaml(sim_inputfile)
    except FileNotFoundError:
        print(f'Failed to load {sim_inputfile}\n', usage)
        raise

    try:
        opts, _ = getopt.getopt(
            argv[2:],
            "h",
            []
        )
    except getopt.GetoptError:
        print(usage)
        raise

    for opt, _ in opts:
        if opt == '-h':
            print(usage)
            sys.exit()

    return calc_params, sim_params
