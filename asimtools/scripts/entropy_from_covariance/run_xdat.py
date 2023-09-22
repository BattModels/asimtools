#!/usr/bin/env python
'''
Runs a user defined lammps script or template 

Author: mkphuthi@github.com
'''
from typing import Dict, Optional, TypeVar
from pathlib import Path
import subprocess
import numpy as np
from ase import Atoms
from ase.io import write
from ase.geometry.geometry import wrap_positions
from asimtools.utils import (
    get_atoms,
    get_images,
)

Atoms = TypeVar('Atoms')

def read_xdat_xyz(xyz_path: str):
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

def write_xdat_xyz(atoms: Atoms, fname: str = 'xyz.ideal'):
    txt = ''
    # atoms.center(about=0)
    # print('pos1:', atoms.positions)
    cell = atoms.get_cell()
    for row in cell:
        txt += f'    {row[0]:.16f}\t{row[1]:.16f}\t{row[2]:.16f}\n'

    txt += f'{len(atoms)} atoms from POSCAR\n'
    symbols_count = {symbol: 1 for symbol in atoms.symbols}
    normalized_cell = (cell.transpose() / np.linalg.norm(cell, axis=1)).transpose()
    # print('normalized_cell:', normalized_cell)
    coords = atoms.get_scaled_positions(wrap=False)
    # print('coords:', coords)
    wrapped_coords = wrap_positions(
        coords,
        normalized_cell,
        pbc=True,
        center=(0.0, 0.0, 0.0),
        pretty_translation=False,
        eps=1e-7,
    )
    # print('wrapped:', wrapped_coords)
    for i, coord in enumerate(wrapped_coords):
        symbol = atoms[i].symbol
        number = atoms[i].number
        # print(coord)
        txt += f'{coord[0]:.10f} {coord[1]:.10f} {coord[2]:.10f} '
        txt += f'{number} {i+1} {symbols_count[symbol]}\n'
        symbols_count[symbol] += 1

    # print(txt)
    with open(fname, 'w', encoding='utf-8') as f:
        f.write(txt)

def write_Mass(atoms: Atoms):
    txt = f'{len(set(atoms.symbols))}\n'
    for number in set(atoms.numbers):
        txt += f'{number} '
    txt += '\n'
    for mass in set(atoms.get_masses()):
        txt += f'{mass} '

    with open('Mass', 'w', encoding='utf-8') as f:
        f.write(txt)


def write_symmetry_ops(atoms: Atoms):
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
    unit_cell_image: Optional[Dict] = None,
    ideal_image: Optional[Dict] = None,
    md_images: Optional[Dict] = None,
    temperature: Optional[float] = None,
) -> Dict:
    '''
    Script does xyz specifically
    '''

    if unit_cell_image is not None:
        unit_cell_atoms = get_atoms(**unit_cell_image)
        write_symmetry_ops(unit_cell_atoms)

        write_Mass(unit_cell_atoms)

    if ideal_image is not None:
        ideal_atoms = get_atoms(**ideal_image)
        write_xdat_xyz(ideal_atoms)

    if md_images is not None:
        images = get_images(**md_images)
        write('XDATCAR', images, format='vasp-xdatcar')

    if temperature is not None:
        with open('Trun', 'w', encoding='utf-8') as f:
            f.write(f'{temperature:.1f}')

    return {}

def run_xdat(
    unit_cell_image: Optional[Dict] = None,
    ideal_image: Optional[Dict] = None,
    md_images: Optional[Dict] = None,
    temperature: Optional[float] = None,
    xyz_ideal: str = 'xyz.ideal',
    xdatcar: str = 'XDATCAR',
    sym_file: str = 'symmetry_operation.dat',
    system: str = 'mod'
) -> Dict:
    ''' 
    Runs the entropy calculation using xdat package
    '''

    gen_xdat_inputs(unit_cell_image, ideal_image, md_images, temperature)

    assert Path(xyz_ideal).exists(), 'Provide existing xyz file'
    assert Path(xdatcar).exists(), 'Provide existing xdatcar file'
    assert Path(sym_file).exists(), 'Provide existing symmetry operation file'

    if system == 'mod':
        command = ['xdat-ent-rhoALL_quantum_mod']
    elif system == 'alloy':
        command = ['xdat-ent-rhoALL_quantum_alloy']
    elif system == 'npt':
        command = ['xdat-ent-rhoALL_quantum_npt']
    else:
        assert False, ['Choose system as one of mod, alloy or npt']

    command += ['-f', xyz_ideal]
    command += ['-s', sym_file]
    command += [xdatcar]

    completed_process = subprocess.run(
            command, check=False, capture_output=True, text=True,
        )

    with open('xdat_stdout.txt', 'w', encoding='utf-8') as f:
        f.write(completed_process.stdout)

    if completed_process.returncode != 0:
        err_txt = f"Failed to run {' '.join(command)}\n"
        err_txt += 'See lmp.stderr.txt for details.'
        print(err_txt)
        with open('xdat_stderr.txt', 'w', encoding='utf-8') as f:
            f.write(completed_process.stderr)
        completed_process.check_returncode()
        return None

    lines = completed_process.stdout.split('\n')

    for line in lines:
        if line.startswith('NPAIR'):
            npair = int(line.split('=')[-1])
        elif line.startswith('Lambda'):
            lambda_val = float(line.split('=')[-1])
        elif line.startswith('T='):
            vals = line.split()
            temperature = float(vals[0].split('=')[-1])
            ln_det = float(vals[1].split('=')[-1])
            s_classical = float(vals[2].split('=')[-1])
            s_quantum = float(vals[3].split('=')[-1])

    results = {
        'npair': npair,
        'lambda': lambda_val,
        'temperature': temperature,
        'ln_det': ln_det,
        's_classical': s_classical,
        's_quantum': s_quantum,
    }

    return results
