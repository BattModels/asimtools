'''
Utilities and helper functions for reading and writing data using set
standards
'''
from typing import (
    TypeVar, Iterable, Tuple, Union, List, Dict, Sequence, Optional
)
from copy import deepcopy
from glob import glob
import os
import logging
import subprocess
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
from ase.io import read
from ase.parallel import paropen
import ase.db
import ase.build
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

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

    if output is None:
        return {}
    return output

def write_yaml(yaml_path: str, yaml_Dict: Dict) -> None:
    """Write a dictionary to a yaml file

    :param yaml_path: Path to write yaml to
    :type yaml_path: str
    :param yaml_Dict: Dictionary to write
    :type yaml_Dict: Dict
    """ 
    # Use paropen so that only the master process is updating outputs
    with paropen(yaml_path, 'w', encoding='utf-8') as f:
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
    header: str = '',
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

    csv_data = pd.DataFrame(data, columns=columns, **kwargs)
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
    name = '__'.join(final_substrs) + '__'
    return name

def get_atoms(
    image_file: Optional[str] = None,
    interface: str = 'ase',
    builder: Optional[str] = 'bulk',
    atoms: Optional[Atoms] = None,
    repeat: Optional[Tuple[int, int, int]] = None,
    rattle_stdev: Optional[float] = None,
    mp_id: Optional[str] = None,
    user_api_key: Optional[str] = None,
    return_type: str = 'ase',
    **kwargs
) -> Union[Atoms, Structure]:
    """Return an ASE Atoms or pymatgen Structure object based on specified 
    config. This is the recommended way to load atoms objects for asimmodules.  

    :param image_file: Path to an ASE-readable image file, defaults to None
    :type image_file: str, optional
    :param interface: Whether to use "ase" or "pymatgen" to create the atoms 
        object
    :type interface: str, optional
    :param builder: Builder to use from :mod:`ase.build`, defaults to 'bulk'. \
        Any extra keyword arguments specified are passed to the build
    :type builder: str, optional
    :param atoms: :class:`ase.Atoms` object, defaults to None
    :type atoms: Atoms, optional
    :param repeat: Number of times to repeat input to create supercell, \
        defaults to None
    :type repeat: Tuple[int, int, int], optional
    :param rattle_stdev: `stdev` to be provided to :meth:`ase.Atoms.rattle`, \
        or as `distance` to :meth:`pymatgen.core.structure.Structure.perturb`
        defaults to None
    :type rattle_stdev: float, optional
    :param mp_id: Material Project ID to fetch structure from, defaults to None
    :type mp_id: str, optional
    :param user_api_key: Material Project API key, must be provided to get
        structures from Materials Project, defaults to None
    :type user_api_key: str, optional
    :param return_type: When set to `ase` returns a :class:`ase.Atoms` object,
        when set to `pymatgen` returns a 
        :class:`pymatgen.core.structure.Structure` object, defaults to 'ase'
    :return: One :class:`ase.Atoms` or 
        :class:`pymatgen.core.structure.Structure` instance
    :rtype: Union[Atoms, Structure]

    There are three options one could use to specify and image or an atoms
    objects:

    #. image_file + \*\*kwargs
    #. builder + \*\*kwargs.
    #. atoms

    Examples
    --------

    Some examples using builders from ASE. All ``**kwargs`` are passed to
    :func:`ase.build`:

    >>> get_atoms(builder='molecule', name='H2O')
    Atoms(symbols='OH2', pbc=False)
    >>> get_atoms(builder='bulk', name='Cu')
    Atoms(symbols='Cu', pbc=True, cell=[[0.0, 1.805, 1.805], [1.805, 0.0, 1.805], [1.805, 1.805, 0.0]])
    >>> get_atoms(builder='bulk', name='Ar', crystalstructure='fcc', a=3.4, cubic=True)
    Atoms(symbols='Ar4', pbc=True, cell=[3.4, 3.4, 3.4])
    >>> get_atoms(builder='fcc100', symbol='Fe', vacuum=8, size=[4,4, 5])
    Atoms(symbols='Cu80', pbc=[True, True, False], cell=[10.210621920333747, 10.210621920333747, 23.22], tags=...)

    Some examples for reading an image from a file using :func:`ase.io.read`
    are given below. All ``**kwargs`` are passed to :func:`ase.io.read`

    >>> h2o = get_atoms(builder='molecule', name='H2O')
    >>> h2o.write('h2o.cif')
    >>> get_atoms(image_file='h2o.cif')
    Atoms(symbols='OH2', pbc=False)
    >>> get_atoms(image_file='h2o.cif', format='cif')
    Atoms(symbols='OH2', pbc=False)
    >>> from ase.io import write
    >>> molecules = [get_atoms(builder='molecule', name='H2O'), get_atoms(builder='molecule', name='H2')]
    >>> write('molecules.xyz', molecules, format='extxyz')
    >>> get_atoms(image_file='molecules.xyz', index=0) # Pick out one structure using indexing
    Atoms(symbols='OH2', pbc=False)

    You can also make supercells and rattle the atoms

    >>> li_bulk = get_atoms(name='Li')
    >>> li_bulk.write('POSCAR', format='vasp')
    >>> get_atoms(image_file='POSCAR', repeat=[3,3,3])
    Atoms(symbols='Li27', pbc=True, cell=[[-5.235, 5.235, 5.235], [5.235, -5.235, 5.235], [5.235, 5.235, -5.235]])
    >>> get_atoms(builder='bulk', name='Li', repeat=[2,2,2], rattle_stdev=0.01)
    Atoms(symbols='Li8', pbc=True, cell=[[-3.49, 3.49, 3.49], [3.49, -3.49, 3.49], [3.49, 3.49, -3.49]])

    Mostly for internal use and use in asimmodules, one can specify atoms
    directly

    >>> li_bulk = get_atoms(name='Li')
    >>> get_atoms(atoms=li_bulk)
    Atoms(symbols='Li', pbc=True, cell=[[-1.745, 1.745, 1.745], [1.745, -1.745, 1.745], [1.745, 1.745, -1.745]])

    In an asimmodule, the ``image`` argument is always given as a dictionary,
    you therefore have to expand it before passing it to ``get_atoms``

    >>> image = {'name': 'Pt'}
    >>> get_atoms(**image)
    Atoms(symbols='Pt', pbc=True, cell=[[0.0, 1.96, 1.96], [1.96, 0.0, 1.96], [1.96, 1.96, 0.0]])

    To get a pymatgen structure object, set the `return_type` to "pymatgen"

    >>> image = {'name': 'Pt', 'return_type': 'pymatgen'}
    >>> get_atoms(**image)
    Structure Summary
    Lattice
    ...

    To download a structure from the Materials Project, you need to provide
    an API key, the MP ID of the structure and set the interface to 'pymatgen'
    You can also specify whether you want the primitive(default) or 
    conventional unit cell as a keyword argument

    >>> {'mp_id': 'mp-14', 'interface': 'pymatgen', 'user_api_key': "USER_API_KEY", 'conventional_unit_cell': True},
    >>> get_atoms(**image)
    Structure Summary
    Lattice
    abc : 2.7718585822512662 2.7718585822512662 2.7718585822512662
    ...
    """
    if interface == 'ase':
        assert image_file is not None or \
            len(kwargs) > 0 or \
            atoms is not None, \
            "Specify atoms, image_file or use ase.build tools if using \
                ASE interface"
        if image_file is not None:
            assert Path(image_file).exists(), \
                f'The image_file {image_file} does not exist from {os.getcwd()}'
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
    
    if interface == 'pymatgen':
        if mp_id is not None:
            from pymatgen.ext.matproj import MPRester
            with MPRester(user_api_key) as mpr:
                struct = mpr.get_structure_by_material_id(mp_id, **kwargs)   
        elif builder is not None:
            builder_func = getattr(Structure, builder)
            try:
                struct = builder_func(**kwargs)
            except ValueError:
                err_txt = 'Make sure you choose a valid builder and '
                err_txt += 'appropriate arguments matching the functions '
                err_txt += 'from pymatgen.core.structure \n'
                err_txt += 'See https://pymatgen.org/pymatgen.core.html'
                err_txt += '#module-pymatgen.core.structure'
                print(err_txt)
                raise

    if repeat is not None:
        if interface == 'ase':
            atoms = atoms.repeat(repeat)
        elif interface == 'pymatgen':
            struct = struct.make_supercell(repeat)

    if rattle_stdev is not None and interface == 'ase':
        atoms.rattle(stdev=rattle_stdev)
    elif rattle_stdev is not None and interface == 'pymatgen':
        struct.perturb(distance=rattle_stdev, min_distance=0)

    if return_type == 'ase' and interface == 'ase':
        return atoms
    elif return_type == 'pymatgen' and interface == 'ase':
        return AseAtomsAdaptor.get_structure(atoms)
    elif return_type == 'pymatgen' and interface == 'pymatgen':
        return struct
    elif return_type == 'ase' and interface == 'pymatgen':
        return AseAtomsAdaptor.get_atoms(struct, msonable=False)

def parse_slice(value: str) -> slice:
    """Parses a :func:`slice` from string, like `start:stop:step`.

    :param value: Slice string
    :type value: str
    :return: slice object
    :rtype: slice
    """
    indices = str(value).split(':')
    parts = []
    for index in indices:
        if index == '':
            parts.append(False)
        else:
            parts.append(index)
    return slice(*[int(p) if p else None for p in parts])

def get_images(
    image_file: str = None,
    pattern: str = None,
    patterns: List[str] = None,
    images: Iterable[Atoms] = None,
    index: Union[str, int] = ':',
    skip_failed: bool = False,
    **kwargs
) -> List[Atoms]:
    """Return a list of atoms objects based on the input arguments. Options \
        to specify are:
        #. image_file
        #. pattern
        #. patterns
        #. images

    :param image_file: Path to ASE-readable file with one or more images, \
        defaults to None
    :type image_file: str, optional
    :param pattern: String pattern of paths from which to search for images, \
        defaults to None. This only gets the last image from each file as in \
        :func:`ase.io.read` if an index is not specified.
    :type pattern: str, optional
    :param patterns: Sequence of string patterns/paths from which to search \
        for images, defaults to None. This only gets one image from each file \
        as in :func:`ase.io.read` without specifying an index
    :type pattern: str, optional
    :param images: A list of atoms objects, defaults to None
    :type images: Iterable[Atoms], optional
    :param index: Index to specify when using :func:`ase.io.read`, \ defaults
        to ':'
    :type index: Union[str, int], optional
    :param skip_failed: Whether to raise an IO error if it fails to read any of
        the specified images or ignore errors, defaults to False
    :type skip_failed: bool, optional
    :raises IOError: Failed to read one of the specified images
    :return: List of :class:`ase.Atoms` for all images found
    :rtype: List[Atoms]

    There are three options one could use to specify and image or an atoms
    objects:

    #. image_file + \*\*kwargs for specifying one image file
    #. pattern + \*\*kwargs for specifying multiple image files with a wildcard
       character
    #. patterns + \*\*kwargs for specifying a list of patterns to match,
       captures the above two cases
    #. images

    Examples
    --------

    Some examples for reading images selectively from an image_file. All
    ``**kwargs`` are passed to :func:`ase.io.read`:

    >>> from asimtools.utils import get_atoms
    >>> molecules = []
    >>> molecules.append(get_atoms(builder='molecule', name='H2O'))
    >>> molecules.append(get_atoms(builder='molecule', name='H2'))
    >>> molecules.append(get_atoms(builder='molecule', name='N2'))
    >>> write('molecules.xyz', molecules, format='extxyz')
    >>> get_images(image_file='molecules.xyz')
    [Atoms(symbols='OH2', pbc=False), Atoms(symbols='H2', pbc=False), Atoms(symbols='N2', pbc=False)]
    >>> get_images(image_file='molecules.xyz', index=':2')
    [Atoms(symbols='OH2', pbc=False), Atoms(symbols='H2', pbc=False)]

    You can also use a wildcard (\*) by specifying the pattern argument. Notice
    that the files don't have to be the same format if ASE can guess all the
    file formats, otherwise you can specify the format argument which should
    apply to all the images.

    >>> cu = get_atoms(name='Cu')
    >>> cu.write('bulk_cu.cfg')
    >>> fe = get_atoms(name='Fe')
    >>> fe.write('bulk_fe.cif')
    >>> pt = get_atoms(name='Pt')
    >>> pt.write('bulk_pt.cfg')
    >>> get_images(pattern='bulk*')
    [Atoms(symbols='Cu', pbc=True, cell=[[0.0, 1.805, 1.805], [1.805, 0.0, 1.805], [1.805, 1.805, 0.0]], masses=..., momenta=...), Atoms(symbols='Fe', pbc=True, cell=[[2.48549, 0.0, 0.0], [-0.8284876429214074, 2.3433456351179887, 0.0], [-0.8284876429214074, -1.171653675382785, 2.0294079014797743]], spacegroup_kinds=...), Atoms(symbols='Pt', pbc=True, cell=[[0.0, 1.96, 1.96], [1.96, 0.0, 1.96], [1.96, 1.96, 0.0]], masses=..., momenta=...)]
    Atoms(symbols='OH2', pbc=False)
    >>> get_images(pattern='bulk*.cfg', format='cfg')
    [Atoms(symbols='Cu', pbc=True, cell=[[0.0, 1.805, 1.805], [1.805, 0.0, 1.805], [1.805, 1.805, 0.0]], masses=..., momenta=...), Atoms(symbols='Pt', pbc=True, cell=[[0.0, 1.96, 1.96], [1.96, 0.0, 1.96], [1.96, 1.96, 0.0]], masses=..., momenta=...)]
    
    You can also specify multiple patterns

    >>> get_images(patterns=['bulk*.cfg', 'bulk\*.cif'])
    [Atoms(symbols='Cu', pbc=True, cell=[[0.0, 1.805, 1.805], [1.805, 0.0, 1.805], [1.805, 1.805, 0.0]], masses=..., momenta=...), Atoms(symbols='Pt', pbc=True, cell=[[0.0, 1.96, 1.96], [1.96, 0.0, 1.96], [1.96, 1.96, 0.0]], masses=..., momenta=...), Atoms(symbols='Fe', pbc=True, cell=[[2.48549, 0.0, 0.0], [-0.8284876429214074, 2.3433456351179887, 0.0], [-0.8284876429214074, -1.171653675382785, 2.0294079014797743]], spacegroup_kinds=...)]
    
    Or you can directly pass a list of Atoms, mostly for internal use

    >>> get_images(images=molecules)
    [Atoms(symbols='OH2', pbc=False), Atoms(symbols='H2', pbc=False), Atoms(symbols='N2', pbc=False)]

    In an asimmodule, the ``images`` argument is always given as a dictionary,
    you therefore have to expand it before passing it to ``get_images``

    >>> images = {'pattern': 'bulk*'}
    >>> get_images(**images)
    [Atoms(symbols='Cu', pbc=True, cell=[[0.0, 1.805, 1.805], [1.805, 0.0, 1.805], [1.805, 1.805, 0.0]], masses=..., momenta=...), Atoms(symbols='Fe', pbc=True, cell=[[2.48549, 0.0, 0.0], [-0.8284876429214074, 2.3433456351179887, 0.0], [-0.8284876429214074, -1.171653675382785, 2.0294079014797743]], spacegroup_kinds=...), Atoms(symbols='Pt', pbc=True, cell=[[0.0, 1.96, 1.96], [1.96, 0.0, 1.96], [1.96, 1.96, 0.0]], masses=..., momenta=...)]
    """
    assert (image_file is not None) or \
        (pattern is not None) or \
        (patterns is not None) or \
        (images is not None), \
        "Please specify file, pattern or iterable"

    if image_file is not None:
        image_file = Path(image_file).resolve()
        images = read(image_file, index=index, **kwargs)
        if not isinstance(images, list):
            images = [images]
    elif pattern is not None:
        image_files = sorted(glob(pattern))
        assert len(image_files) > 0, \
            f'No images matching pattern "{pattern}" from "{os.getcwd()}"'

        images = []
        for image_file in image_files:
            image_file = Path(image_file).resolve()
            try:
                new_images = read(
                    image_file,
                    index=index,
                    **kwargs
                )
            except Exception as exc:
                if not skip_failed:
                    raise IOError(
                        f"Failed to read {image_file} from {os.getcwd()}"
                    ) from exc
                else:
                    new_images = []

            # Output of read can either be list of atoms or Atoms, depending on index
            if not isinstance(new_images, list):
                new_images = [new_images]
            images += new_images

    elif patterns is not None:
        images = []
        for pattern in patterns:
            image_files = sorted(glob(pattern))
            assert len(image_files) > 0, \
                f'Don\'t include pattern "{pattern}" if no files match'
            images += get_images(
                pattern=pattern,
                index=index,
                skip_failed=skip_failed,
                **kwargs
            )
    elif images is not None:
        images = images[parse_slice(str(index))]
    else:
        images = []

    if not skip_failed:
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
    then ``.asimtools`` dotfile. If an env_input.yaml is present in the 
    work direcory, it takes precedence.

    :return: env_input dictionary or empty dictionary if not found
    :rtype: Dict
    """
    if Path('env_input.yaml').exists():
        env_input_file = Path('env_input.yaml').resolve()
    else:
        env_input_file = os.getenv('ASIMTOOLS_ENV_INPUT', None)

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
                except FileNotFoundError:
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
                except FileNotFoundError:
                    print('WARNING: No calc_input yaml found')
    else:
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

def change_dict_value(
    d: Dict,
    new_value,
    key_sequence: Sequence,
    return_copy: bool  = True
) -> Dict:
    """Changes a value in the specified dictionary given by following the 
    key sequence

    :param d: dictionary to be changed
    :type d: Dict
    :param new_value: The new value that will replace the old one
    :type new_value: _type_
    :param key_sequence: List of keys in the order in which they access the dictionary key
    :type key_sequence: Sequence
    :param return_copy: Whether to return a copy only or to modify the dictionary in-place as well, defaults to True
    :type return_copy: bool, optional
    :return: The changed dictionary
    :rtype: Dict
    """
    if return_copy:
        d = deepcopy(d)
    if len(key_sequence) == 1:
        d[key_sequence[0]] = new_value
        return d
    else:
        new_d = change_dict_value(
            d[key_sequence[0]],
            new_value,
            key_sequence[1:],
            return_copy=return_copy
        )
        d[key_sequence[0]] = new_d
        return d

def change_dict_values(
    d: Dict,
    new_values: Sequence,
    key_sequences: Sequence,
    return_copy: bool = True
) -> Dict:
    """Changes values in the specified dictionary given by following the 
    key sequences. Key-value pairs are set in the given order

    :param d: dictionary to be changed
    :type d: Dict
    :param new_values: The new values that will replace the old one
    :type new_values: Sequence
    :param key_sequence: List of list of keys in the order in which they access the dictionary key
    :type key_sequence: Sequence
    :param return_copy: Whether to return a copy only or to modify the dictionary in-place as well, defaults to True
    :type return_copy: bool, optional
    :return: The changed dictionary
    :rtype: Dict
    """
    for key_sequence, new_value in zip(key_sequences, new_values):
        d = change_dict_value(d, new_value, key_sequence, return_copy)

    return d

def get_logger(
    logfile='job.log',
    fstr='%(asctime)s |%(module)s|%(funcName)s| %(levelname)s: %(message)s',
    level='debug',
):
    ''' Get the logger '''

    level_dict = {
        'debug': logging.DEBUG,
        'info': logging.INFO,
        'warning': logging.WARNING,
    }
    assert level in level_dict, f'level must be one of {level_dict.keys()}'

    logging.basicConfig(
        filename=logfile,
        level=level_dict[level],
        format=fstr,
    )
    return logging.getLogger(__name__)


def get_str_btn(
    s: Union[str,os.PathLike],
    s1: str,
    s2: str,
    occurence: Optional[int] = 0,
    start_index: Optional[int] = 0,
):
    """Returns the substring between strings s1 and s2 from s

    :param s: string/path from which to extract substring
    :type s: str
    :param s1: substring before the desired substring, None starts from 
        the beginning of s
    :type s1: str
    :param s2: substring after the desired substring, None ends at the 
        end of the string
    :type s2: str
    :param occurence: which occurence to return e.g. occurence=1 returns the 
        second time a substring is encountered, defaults to 0
    :type occurence: int, optional
    :param start_index: first index from which to search for substring, 
        defaults to 0
    :type int: int, optional
    :return: substring
    :rtype: _type_
    """

    j = 0
    stop_index = len(s) + 1

    s = s[start_index:stop_index]
    while occurence - j >= 0:
        if s1 is not None:
            i1 = s.index(s1) + len(s1)
        else:
            i1 = 0
        if s2 is not None:
            i2 = s[i1:].index(s2) + i1
        else:
            i2 = len(s)

        if occurence - j == 0:
            return s[i1:i2]

        s = s[i1:]
        j += 1

def find_nth(haystack: str, needle: str, n: int) -> int:
    ''' Return index of nth occurence of substring in string '''
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start

def get_nth_label(
    s: str,
    n: int = 1,
):
    ''' Return nth label in a string potentially containing multiple labels,
    indexing starts from 0 '''
    start = find_nth(s, '__', n=(n*2+1))
    return get_str_btn(s, '__', '__', start_index=start)
