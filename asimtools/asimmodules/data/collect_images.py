'''
Collects images into a single file

author: mkphuthi@github.com
'''

from typing import Dict, Optional, Sequence
import numpy as np
from ase.io import write
from asimtools.utils import (
    get_atoms, get_images, new_db
)


def collect_images(
    images: Dict,
    out_format: str = 'extxyz',
    fnames: Sequence[str] = ['output_images.xyz'],
    splits: Optional[Sequence[float]] = (1,),
    shuffle: bool = True,
    rename_keys: Optional[Dict] = None,
    energy_per_atom_limits: Optional[Sequence[float]] = None,
    force_max: Optional[float] = None,
    stress_limits: Optional[Sequence[float]] = None,
) -> Dict:
    """Collects images into one file/database and can split them into 
    multiple files/databases for ML tasks

    :param images: Images specification, see :func:`asimtools.utils.get_images`
    :type images: Dict
    :param out_format: output file format, defaults to 'extxyz'
    :type out_format: str, optional
    :param fnames: file name without extension, defaults to 'output_images.xyz'
    :type fnames: str, optional
    :param splits: Ratios to split data into, defaults to None
    :type splits: Optional[Sequence[float]], optional
    :param shuffle: shuffle images before splitting, defaults to True
    :type shuffle: bool, optional
    :param rename_keys: keys to rename on writing to the output file, defaults to None
    :type rename_keys: Optional[Dict], optional
    :param energy_per_atom_limits: energy limits for filtering images, defaults to None
    :type energy_per_atom_limits: Optional[Sequence[float]], optional
    :param force_max: forces maximimum for filtering images, defaults to None
    :type force_max: Optional[float], optional
    :param stress_limits: stress limits for filtering images, defaults to None
    :type stress_limits: Optional[Sequence[float]], optional
    :return: results
    :rtype: Dict
    
    """
    write_kwargs = {}
    if fnames == (1):
        fnames = [f'{fnames}-{i:03d}' for i in range(len(splits))]
    assert len(splits) == len(fnames), \
        'Number of splits must match number of fnames'
    images = get_images(**images)

    selected_atoms = []
    nonselected_atoms = []
    for atoms in images:
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        stress = atoms.get_stress()

        select = True
        if energy_per_atom_limits is not None:
            if (energy < energy_per_atom_limits[0]):
                if (energy > energy_per_atom_limits[1]):
                    select = False
        if force_max is not None:
            max_force = np.max(np.linalg.norm(forces, axis=1))
            if (max_force > force_max):
                select = False
        if stress_limits is not None:
            max_stress = np.max(stress)
            min_stress = np.min(stress)
            if (min_stress < stress_limits[0]):
                if (max_stress > stress_limits[1]):
                    select = False
        print(select)

        if select:
            if rename_keys is not None and out_format == 'extxyz':
                write_kwargs['write_info'] = True
                if 'energy' in rename_keys:
                    atoms.info[rename_keys['energy']] = energy
                if 'forces' in rename_keys:
                    atoms.arrays[rename_keys['forces']] = forces
                if 'stress' in rename_keys:
                    atoms.info[rename_keys['stress']] = stress
            selected_atoms.append(atoms)
        else:
            nonselected_atoms.append(atoms)

    write('nonselected_images.xyz', nonselected_atoms, format='extxyz')
    if shuffle:
        np.random.shuffle(selected_atoms)
    datasets = []
    start_index = 0
    index_ranges = []
    for i, split in enumerate(splits):
        index_range = int(split / np.sum(splits) * len(selected_atoms))
        index_ranges.append(index_range)
        dataset = selected_atoms[start_index:start_index + index_range]
        if format == 'db':
            with new_db(f'{fnames[i]}.db') as db:
                db.write(dataset)
        else:
            write(f'{fnames[i]}', dataset, format=out_format, **write_kwargs)
        start_index += index_range

    return {"split": index_ranges, 'num_rejected': len(nonselected_atoms)}
