'''
Utilities and helper functions for reading and writing data using set
standards
'''
from typing import (
    TypeVar, Iterable, Tuple, Union, List, Dict, Sequence, Optional
)
from glob import glob
import os
import subprocess
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
from ase.io import read
import ase.db
import ase.build

Atoms = TypeVar('Atoms')

def read_yaml(yaml_path: str) -> Dict:
    """Read a yaml file

    :param yaml_path: Path to yaml file
    :type yaml_path: str
    :return: Dictionary
    :rtype: Dict
    """
    with open(yaml_path, 'r', encoding='utf-8') as f:
        output = yaml.safe_load(f)
    return output

def write_yaml(yaml_path: str, yaml_Dict: Dict) -> None:
    """Write a dictionary to a yaml file

    :param yaml_path: Path to write yaml to
    :type yaml_path: str
    :param yaml_Dict: Dictionary to write
    :type yaml_Dict: Dict
    """ 
    with open(yaml_path, 'w', encoding='utf-8') as f:
        yaml.dump(yaml_Dict, f)

def get_axis_lims(x: Sequence, y: Sequence, padding: float=0.1):
    """Get an estimate of good limits for a plot axis"""
    data_min = np.min([x, y])
    data_max = np.max([x, y])
    diff = data_max - data_min
    lims = [data_min - padding * diff, data_max + padding * diff]
    return lims

def write_csv_from_dict(
    fname: str,
    data: Dict,
    columns: Optional[Sequence] = None,
    header: Optional[str] = None,
    **kwargs
) -> pd.DataFrame:
    """Write a tabulated csv to a file from a dictionary. The keys will
    become the columns and the values will be the columns, All columns used 
    must have the same length. kwargs are passed to 
    :func:`pandas.Dataframe.to_csv()`

    :param fname: File name to be return to
    :type fname: str
    :param data: Dictionary to write
    :type data: Dict
    :param columns: Keys to use as columns of csv, defaults to None
    :type columns: Iterable, optional
    :return: Pandas dataframe of results
    :rtype: pd.DataFrame
    """
    if columns is None:
        columns = list(data.keys())

    csv_data = pd.DataFrame(data, columns=columns)
    with open(fname, 'w', encoding='utf-8') as f:
        f.write('#' + header + '\n')

    csv_data.to_csv(fname, index=False, header=True, mode='a')
    return csv_data

def strip_symbols(substr: str) -> str:
    """Helper function to format filenames using standard

    :param substr: substring
    :type substr: str
    :return: stripped string
    :rtype: str
    """
    bad_symbols = ['_', '-', '.', ' ']
    if len(substr) == 0:
        return ''
    if substr[0] in bad_symbols or substr[-1] in bad_symbols:
        for bad_symbol in bad_symbols:
            substr = substr.strip(bad_symbol)
        return strip_symbols(substr)
    return substr

def join_names(substrs: Sequence[str]) -> str:
    """Join multiple strings using a hyphen neatly

    :param substrs: Sequence of strings to merge
    :type substrs: Sequence[str]
    :return: Joined strings
    :rtype: str
    """
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
    """Return an atoms object based on specified config. This is the 
    recommended way to load atoms objects. There are three options to specify:
    #. image_file.
    #. builder, kwargs.
    #. atoms.

    :param image_file: Path to an ASE-readable image file, defaults to None
    :type image_file: str, optional
    :param builder: Builder to use from :mod:`ase.build`, defaults to 'bulk'. \
        Any extra keyword arguments specified are passed to the build
    :type builder: str, optional
    :param atoms: :class:`ase.Atoms` object, defaults to None
    :type atoms: Atoms, optional
    :param repeat: Number of times to repeat input to create supercell, \
        defaults to None
    :type repeat: Tuple[int, int, int], optional
    :param rattle_stdev: stdev to be provided to :meth:`ase.Atoms.rattle`, \
        defaults to None
    :type rattle_stdev: float, optional
    :return: One :class:`ase.Atoms` instance
    :rtype: Atoms
    """
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

def parse_slice(value: str) -> slice:
    """Parses a :func:`slice` from string, like `start:stop:step`.

    :param value: Slice string
    :type value: str
    :return: slice object
    :rtype: slice
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
    patterns: List[str] = None,
    images: Iterable[Atoms] = None,
    index: Union[str, int] = ':',
    **kwargs
) -> List[Atoms]:
    """Return a list of atoms objects based on the input arguments. Options \
        to specify are:
        #. image_file
        #. pattern
        #. images

    :param image_file: Path to ASE-readable file with one or more images, \
        defaults to None
    :type image_file: str, optional
    :param pattern: String pattern of paths from which to search for images, \
        defaults to None. This only gets one image from each file as in \
        :func:`ase.io.read` without specifying an index
    :type pattern: str, optional
    :param patterns: Sequence of string patterns/paths from which to search \
        for images, defaults to None. This only gets one image from each file \
        as in :func:`ase.io.read` without specifying an index
    :type pattern: str, optional
    :param images: A list of atoms objects, defaults to None
    :type images: Iterable[Atoms], optional
    :param index: Index to specify when using :func:`ase.io.read`, \
        defaults to ':'
    :type index: Union[str, int], optional
    :return: List of :class:`ase.Atoms`
    :rtype: List[Atoms]
    """
    assert (image_file is not None) or \
        (pattern is not None) or \
        (images is not None), \
        "Please specify file, pattern or iterable"

    index = parse_slice(str(index))
    if image_file is not None:
        images = read(image_file, index=index, **kwargs)
    elif pattern is not None:
        image_files = glob(pattern)
        assert len(image_files) > 0, \
            f'No images matching pattern "{pattern}" from "{os.getcwd()}"'
        images = [
            read(image_file, **kwargs) for image_file in image_files
        ]
    elif patterns is not None:
        images = []
        for pattern in patterns:
            image_files = glob(pattern)
            assert len(image_files) > 0, \
                f'No images matching patterns "{pattern}" from "{os.getcwd()}"'
            images += get_images(pattern=pattern, **kwargs)
    elif images is not None:
        images = images[index]
    else:
        images = []

    assert len(images) > 0, 'No images found'

    return images

def new_db(dbname: str):
    """Creates a new ASE database overwriting an old one

    :param dbname: Name for the new database
    :type dbname: str
    :return: Connection to database
    :rtype: :mod:`ase.db` connection
    """
    with open(dbname, 'w', encoding='utf-8') as f:
        f.write("")
    return ase.db.connect(dbname)

def get_env_input() -> Dict:
    """Gets the global env_input from the environment variable, then tries
    then ``.asimtools`` dotfile

    :return: env_input dictionary or empty dictionary if not found
    :rtype: Dict
    """
    env_input_file = os.getenv('ASIMTOOLS_ENV_INPUT')
    if env_input_file is None:
        dotfile = Path('~/.asimtools')
        if dotfile.exists():
            with open(dotfile, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            for line in lines:
                if line.startswith('ENV_INPUT_FILE'):
                    env_input_file = Path(line.split('='))
                    break
            if env_input_file.exists():
                try:
                    env_input = read_yaml(env_input_file)
                    return env_input
                except Exception:
                    print('WARNING: No env_input yaml found')
    else:
        env_input = read_yaml(env_input_file)
        return env_input
    return {}

def get_calc_input():
    """Gets the global calc_input from the environment variable, then tries
    then ``.asimtools`` dotfile

    :return: calc_input dictionary or empty dictionary if not found
    :rtype: Dict
    """
    calc_input_file = os.getenv('ASIMTOOLS_CALC_INPUT')
    if calc_input_file is None:
        dotfile = Path('~/.asimtools')
        if dotfile.exists():
            with open(dotfile, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            for line in lines:
                if line.startswith('CALC_INPUT_FILE'):
                    calc_input_file = Path(line.split('='))
                    break
            if calc_input_file.exists():
                try:
                    calc_input = read_yaml(calc_input_file)
                    return calc_input
                except Exception:
                    print('WARNING: No calc_input yaml found')
    else:
        print('calcdfasdfasdf', calc_input_file)
        calc_input = read_yaml(calc_input_file)
        return calc_input
    return {}

def check_if_slurm_job_is_running(slurm_job_id: Union[str,int]):
    """Checks if the slurm job with specifying job_id is running

    :param slurm_job_id: Job ID as str or int
    :type slurm_job_id: Union[str,int]
    :return: Whether job is running or not
    :rtype: bool
    """
    slurm_job_id = str(slurm_job_id)
    completed_process = subprocess.run(
        ['squeue', '--job', slurm_job_id], 
        check=False,
        capture_output=True,
        text=True,
    )
    stdout = completed_process.stdout
    if str(' '+slurm_job_id) in stdout:
        return True
    else:
        return False
    
