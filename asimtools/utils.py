'''
Utilities and helper functions for reading and writing data using set
standards
'''
from typing import TypeVar, Iterable, Tuple, Union
from glob import glob
import argparse
import yaml
import pandas as pd
from ase.io import read
import ase.build

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
    image_file: str = None,
    builder: str = 'bulk',
    atoms: Atoms = None,
    repeat: Tuple[int, int, int] = None,
    rattle_stdev: float = None,
    **kwargs
) -> Atoms:
    '''
    Helper for returning an Atoms object by ASE specification, image
    file or directly using an atoms object
    '''
    assert image_file is not None or \
        len(kwargs) > 0 or \
        atoms is not None, \
        "Specify atoms, image_file or use ase.build tools"

    if image_file is not None:
        atoms = read(image_file, **kwargs)
    elif atoms is None:
        builder_func = getattr(ase.build, builder)
        try:
            atoms = builder_func(**kwargs)
        except ValueError:
            print('Make sure you choose a valid builder and appropriate\
                arguments matching the functions in ase.build. \n\
                See https://wiki.fysik.dtu.dk/ase/ase/build/build.html')
            raise
    else:
        assert atoms is not None, 'Specify an input structure'

    if repeat is not None:
        atoms = atoms.repeat(repeat)

    if rattle_stdev is not None:
        atoms.rattle(stdev=rattle_stdev)

    return atoms

def parse_slice(value):
    """
    Parses a `slice()` from string, like `start:stop:step`.
    """
    if value:
        parts = value.split(':')
        if len(parts) == 1:
            parts = [None, parts[0]]
    else:
        parts = []
    return slice(*[int(p) if p else None for p in parts])

def get_images(
    image_file: str = None,
    pattern: str = None,
    images: Iterable[Atoms] = None,
    index: Union[str, int] = ':',
    **kwargs
) -> list[Atoms]:
    '''
    Helper for returning an Atoms object by ASE specification, input
    file or directly using an atoms object
    '''
    assert (image_file is not None) or \
        (pattern is not None) or \
        (images is not None), \
        "Please specify file, pattern or iterable"

    index = parse_slice(str(index))
    if image_file is not None:
        images = read(image_file, index=index, **kwargs)
    elif pattern is not None:
        image_files = glob(pattern)
        assert len(image_files) > 0, 'No images matching pattern'
        images = [
            read(image_file, **kwargs) for image_file in image_files
        ]
    elif images is not None:
        images = images[index]
    else:
        images = []

    assert len(images) > 0, 'No images found'

    return images

def parse_command_line(args) -> Tuple[dict, dict]:
    ''' Parse command line input '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'calc_input_file',
        metavar='calculator_configuration_file',
        type=str,
        help='calculator configuration yaml file'
    )
    parser.add_argument(
        'sim_input_file',
        metavar='simulation_configuration_file',
        type=str,
        help='calculator configuration yaml file'
    )
    args = parser.parse_args(args)

    calc_params = read_yaml(args.calc_input_file)
    sim_params = read_yaml(args.sim_input_file)

    return calc_params, sim_params
