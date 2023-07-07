'''
Utilities and helper functions for reading and writing data using set
standards
'''
from typing import TypeVar, Iterable, Tuple, Union, List, Dict
from glob import glob
import os
import subprocess
from pathlib import Path
import yaml
import pandas as pd
from ase.io import read
import ase.db
import ase.build

Atoms = TypeVar('Atoms')

def read_yaml(yaml_path: str) -> Dict:
    ''' Read a yaml file'''
    with open(yaml_path, 'r', encoding='utf-8') as f:
        output = yaml.safe_load(f)
    return output

def write_yaml(yaml_path: str, yaml_Dict: Dict) -> None:
    ''' Write a yaml file '''
    with open(yaml_path, 'w', encoding='utf-8') as f:
        yaml.dump(yaml_Dict, f)

def write_csv_from_dict(fname: str, data: Dict, columns: Iterable = None) -> pd.DataFrame:
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
) -> List[Atoms]:
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
        assert len(image_files) > 0, \
            f'No images matching pattern "{pattern}" in "{os.getcwd()}"'
        images = [
            read(image_file, **kwargs) for image_file in image_files
        ]
    elif images is not None:
        images = images[index]
    else:
        images = []

    assert len(images) > 0, 'No images found'

    return images

def new_db(dbname):
    ''' Creates a new ASE database overwriting any existing one '''
    with open(dbname, 'w', encoding='utf-8') as f:
        f.write("")
    return ase.db.connect(dbname)

def get_env_input():
    ''' Gets the global env_input specified by Environment variable '''
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
    ''' Gets the global calc_input specified by environment variable '''
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

def check_if_slurm_job_is_running(slurm_job_id):
    ''' Checks if the job with the given slurm ID is running '''
    completed_process = subprocess.run(
        ['squeue', '--job', str(slurm_job_id)], 
        check=False,
        capture_output=True,
        text=True,
    )
    stdout = completed_process.stdout
    # I know, I know. Not 100% safe but should work most of the time 
    if str(slurm_job_id) in stdout:
        return True
    else:
        return False
    
