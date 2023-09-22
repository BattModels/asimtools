#!/usr/bin/env python
'''
Describe the script briefly here. If this is a script that runs multiple steps,
describe it here using reStructuredText to generate autodocs

Cite the papers where the method/script was first introduced here as well

Author: mkphuthi@github.com
'''
from typing import Dict
import subprocess
import numpy as np
from ase import Atoms
from ase.io import write
from asimtools.utils import (
    get_atoms,
    get_images,
)

def read_xdat_xyz(xyz_path):
    with open(xyz_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    cell = []
    for line in lines[:3]:
        row = [float(val) for val in line.split()]
        cell.append(row)
    cell = np.array(cell)

    natoms = int(lines[3].split()[0])

    # atoms = Atoms(cell=cell, pbc=True)
    scaled_positions = []
    atomic_numbers = []
    for line in lines[4:natoms+4]:
        row = [float(val) for val in line.split()]
        scaled_positions.append(row[:3])
        atomic_numbers.append(int(row[3]))
        # symbol = chemical_symbols[int(row[3])]
        # atom = Atom(symbol=symbol)
        # atom.scaled_position = np.array(scaled_position)
        # atoms.append(atom)
    atoms = Atoms(f'X{natoms}', cell=cell, pbc=True)
    atoms.set_atomic_numbers(atomic_numbers)
    atoms.set_scaled_positions(np.array(scaled_positions))
    return atoms

def write_xdat_xyz(atoms, fname='xyz.ideal'):
    txt = ''
    print(atoms)
    atoms.center()
    cell = atoms.get_cell()
    for row in cell:
        txt += f'    {row[0]:.16f} {row[1]:.16f} {row[2]:.16f}\n'

    txt += f'{len(atoms)} atoms from POSCAR\n'
    symbols_count = {symbol: 1 for symbol in atoms.symbols}
    for i, coord in enumerate(atoms.get_scaled_positions(wrap=True)):
        symbol = atoms[i].symbol
        number = atoms[i].number
        print(coord)
        txt += f'{coord[0]:.10f} {coord[0]:.10f} {coord[0]:.10f} '
        txt += f'{number} {i+1} {symbols_count[symbol]}\n'
        symbols_count[symbol] += 1

    with open(fname, 'w', encoding='utf-8') as f:
        f.write(txt)


def write_Mass(atoms):
    txt = f'{len(set(atoms.symbols))}\n'
    for number in set(atoms.numbers):
        txt += f'{number} '
    txt += '\n'
    for mass in set(atoms.get_masses()):
        txt += f'{mass} '

    with open('Mass', 'w', encoding='utf-8') as f:
        f.write(txt)


def write_symmetry_ops(atoms):
    atoms.write('sym_input.poscar', format='vasp')
    txt = "WARNING: Make sure to use high symmetry unit cell, "
    txt += "not primitive/supercell, e.g. cubic BCC rather than primitive cell"
    print(txt)
    command = "aflow-sym.sh -p sym_input.poscar".split()
    completed_process = subprocess.run(
            command, check=False, capture_output=True, text=True,
        )

    with open('aflow-sym_stdout.txt', 'w', encoding='utf-8') as f:
        f.write(completed_process.stdout)

    if completed_process.returncode != 0:
        err_txt = 'ERROR: See aflow-sym.stderr.txt for details.'
        print(err_txt)
        with open('aflow-sym_stderr.txt', 'w', encoding='utf-8') as f:
            f.write(completed_process.stderr)
        completed_process.check_returncode()
        return None

def gen_xdat_inputs(
    unit_cell_image,
    ideal_image,
    md_images,
    temperature,
) -> Dict:
    '''
    Script does xyz specifically
    '''

    unit_cell_atoms = get_atoms(**unit_cell_image)
    write_symmetry_ops(unit_cell_atoms)

    write_Mass(unit_cell_atoms)

    ideal_atoms = get_atoms(**ideal_image)
    write_xdat_xyz(ideal_atoms)

    images = get_images(**md_images)
    write('XDATCAR', images, format='vasp-xdatcar')

    with open('Trun', 'w', encoding='utf-8') as f:
        f.write(f'{temperature:.1f}')

    results = {}
    return results
