#!/usr/bin/env python
'''
'''

from typing import Dict, Tuple, Sequence
import numpy as np
from phonopy import Phonopy
from phonopy.interface.calculator import read_crystal_structure, write_crystal_structure
from phonopy.structure.atoms import PhonopyAtoms
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms, get_logger

def generate_phonopy_displacements(
    image: Dict,
    supercell: Sequence[int],
    distance: float = 0.01,
    phonopy_save_path: str = 'phonopy_params.yaml'
) -> Dict:

    atoms = get_atoms(**image)
    atoms.write('POSCAR-unitcell', format='vasp')

    unitcell, _ = read_crystal_structure(
        "POSCAR-unitcell",
        interface_mode='vasp'
    )

    phonon = Phonopy(
        unitcell,
        supercell_matrix=[
            [supercell[0], 0, 0],
            [0, supercell[1], 0],
            [0, 0, supercell[2]]
        ],
    )

    phonon.generate_displacements(distance=distance)
    supercells = phonon.supercells_with_displacements

    for i, cell in enumerate(supercells):
        ID = f'{i+1:03d}'
        write_crystal_structure(
            f'supercell-{ID}',
            cell,
            interface_mode='vasp',
            optional_structure_info={'supercell_id': f'{ID}'}
        )

    phonon.save(phonopy_save_path)

    return {}
