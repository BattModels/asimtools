#!/usr/bin/env python
'''
Relaxes the given atomic structure using ASE's built-in structure
optimizers
'''

from typing import Dict, Sequence
from phonopy import Phonopy
from phonopy.interface.calculator import write_crystal_structure
from phonopy.structure.atoms import PhonopyAtoms
from asimtools.utils import get_atoms

def generate_phonopy_displacements(
    image: Dict,
    supercell: Sequence[int],
    distance: float = 0.01,
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
        write_crystal_structure(
            f'POSCAR_disp_{i}',
            cell,
            interface_mode='vasp',
        )

    return {}
