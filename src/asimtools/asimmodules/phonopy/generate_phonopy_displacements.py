#!/usr/bin/env python
'''
'''

from typing import Dict, Tuple, Sequence
import numpy as np
from phonopy import Phonopy
from phonopy.interface.calculator import read_crystal_structure, write_crystal_structure
from phonopy.structure.atoms import PhonopyAtoms
from ase import Atoms
from ase.spacegroup.symmetrize import check_symmetry
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms, get_logger
import spglib

def generate_phonopy_displacements(
    image: Dict,
    supercell: Sequence[int],
    distance: float = 0.01,
    phonopy_save_path: str = 'phonopy_params.yaml',
    refine_cell: bool = False,
    symprec: float = 1e-4,
) -> Dict:
    """Generates displacements for phonopy calculations

    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param supercell: Size of supercell relative to input image
    :type supercell: Sequence[int]
    :param distance: distance to displace atoms, defaults to 0.01
    :type distance: float, optional
    :param phonopy_save_path: Where the phonopy save file is saved, other asimmodule might use it, defaults to 'phonopy_params.yaml'
    :type phonopy_save_path: str, optional
    :return: Results
    :rtype: Dict
    """

    atoms = get_atoms(**image)
    unrefined_sg = check_symmetry(atoms, symprec=symprec).international
    if refine_cell:
        cell = (
            atoms.get_cell().array,
            atoms.get_scaled_positions(),
            atoms.get_atomic_numbers()
        )

        refined_cell = spglib.standardize_cell(cell, symprec=symprec)
        atoms = Atoms(
            cell=refined_cell[0],
            scaled_positions=refined_cell[1],
            numbers=refined_cell[2],
            pbc=True,
        )
        atoms.write('refined_atoms.cif')
        sg = check_symmetry(atoms, symprec=1e-6).international
    else:
        sg = None

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

    results = {
        'unrefined_spacegroup': str(unrefined_sg),
        'refined_spacegroup': str(sg),
    }
    return results
