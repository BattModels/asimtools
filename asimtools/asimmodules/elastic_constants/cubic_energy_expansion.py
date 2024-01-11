#!/usr/bin/env python
'''
Describe the asimmodule briefly here. If this is a asimmodule that runs multiple steps,
describe it here using reStructuredText to generate autodocs

Cite the papers where the method/asimmodule was first introduced here as well

Author: mkphuthi@github.com
'''

from typing import Dict, Sequence, Optional
import logging
import numpy as np
from scipy.optimize import curve_fit
from ase.units import GPa
from asimtools.calculators import load_calc
from asimtools.utils import (
    get_atoms,
)
from asimtools.asimmodules.geometry_optimization.ase_cubic_eos_optimization import (
    ase_cubic_eos_optimization as eos
)

def get_strained_atoms(atoms, strain: str, delta: float):
    """Returns a unit cell with atoms strained according to some useful types

    :param atoms: input atoms
    :type atoms: ase.Atoms      
    :param strain: Type of strain
    :type strain: str
    :param delta: value of strain
    :type delta: float
    :return: atoms with strained cell
    :rtype: ase.Atoms
    """
    strains = {
        'uniform':np.array([
            [1+delta, 0, 0],
            [0, 1+delta, 0],
            [0, 0, 1+delta],
        ]),

        'c44_vol_cons': np.array([
            [4 / (4 - delta**2), 0, 0],
            [0, 1, delta/2],
            [0, delta/2, 1],
        ]),
        'cprime_vol_cons': np.array([
            [1+delta/2, 0, 0],
            [0, 1-delta/2, 0],
            [0, 0, 4 / (4 - delta**2)],
        ]),
        'orth_vol_cons': np.array([
            [1+delta, 0, 0],
            [0, 1-delta, 0],
            [0, 0, 1 + delta**2 / (1 - delta**2)],
        ]),
        'mono_vol_cons': np.array([
            [1, delta/2, 0],
            [delta/2, 1, 0],
            [0, 0, 1+delta**2 / (4 - delta**2)],
        ]),
    }

    if isinstance(strain, str):
        try:
            strain = strains[strain]
        except ValueError:
            print(f' Please choose strain from {strains}')

    cell = atoms.get_cell().copy()
    new_cell = np.array([
        np.squeeze(np.asarray(np.dot(strain, cell[0]))),
        np.squeeze(np.asarray(np.dot(strain, cell[1]))),
        np.squeeze(np.asarray(np.dot(strain, cell[2]))) 
    ])

    strained_atoms = atoms.copy()
    strained_atoms.set_cell(new_cell, scale_atoms=True)
    return strained_atoms

def cubic_energy_expansion(
    calc_id: str,
    image: Dict,
    deltas: Sequence[float] = (-0.01,-0.0075,-0.005,0.00,0.005,0.0075,0.01),
    ase_cubic_eos_args: Optional[Dict] = None,
) -> Dict:
    """Calculates B (Bulk modulus), C11, C12 and C44 elastic constants of 
    a structure with cubic symmetry

    :param calc_id: calc_id specification
    :type calc_id: str
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param deltas: strains to apply in each direction, defaults to (-0.01,-0.0075,-0.005,0.00,0.005,0.0075,0.01)
    :type deltas: Sequence[float], optional
    :param ase_cubic_eos_args: Argumenents to pass to :func:`asimtools.asimmodules.geometry_optimization.ase_cubic_eos_optimization`, defaults to None
    :type ase_cubic_eos_args: Optional[Dict], optional
    :return: Elastic constant results
    :rtype: Dict
    """
    calc = load_calc(calc_id)
    atoms = get_atoms(**image)
    atoms.set_calculator(calc)

    # Start by getting the Bulk modulus and optimized cell from the EOS
    logging.info('Calculating EOS')
    eos_kwargs = {'image': image, 'calc_id': calc_id}
    if ase_cubic_eos_args is not None:
        eos_kwargs.update(ase_cubic_eos_args)
    eos_results = eos(**eos_kwargs)
    B = eos_results['B']
    atoms = get_atoms(image_file='eos_image_output.xyz')
    logging.info('Finished calculating EOS and B, moving to other constants')

    # Sequentially apply relevant strains and extract elastic constants by
    # fitting stress strain curves
    c11min12_ens = []
    c44_ens = []
    for delta in deltas:
        logging.info('Calculating for C44 strain, delta = %s', delta)
        c44_atoms = get_strained_atoms(
            atoms.copy(), 'mono_vol_cons', delta
        )
        calc = load_calc(calc_id)
        c44_atoms.calc = calc
        c44_en = c44_atoms.get_potential_energy()
        c44_ens.append(c44_en)
        c44_atoms.info['strain'] = 'mono_vol_cons'
        c44_atoms.info['delta'] = delta
        c44_atoms.write(f's{delta:.4f}_c44.xyz')

        logging.info('Calculating for C11-C12 strain, delta = %s', delta)
        c11min12_atoms = get_strained_atoms(
            atoms.copy(), 'orth_vol_cons', delta
        )
        calc = load_calc(calc_id)
        c11min12_atoms.calc = calc
        c11min12_en = c11min12_atoms.get_potential_energy()
        c11min12_ens.append(c11min12_en)
        c11min12_atoms.info['strain'] = 'mono_vol_cons'
        c11min12_atoms.info['delta'] = delta
        c11min12_atoms.write(f's{delta:.4f}_c11min12.xyz')

    def f(x, a, b, c):
        ''' Fitting function for free energy expansion to second order'''
        return a*x**2 + b*x + c

    logging.info('Fitting for C44')
    popt, pcov = curve_fit(f, deltas, c44_ens)
    C44 = 2*popt[0] / atoms.get_volume()

    logging.info('Fitting for C11-C12')
    popt, pcov = curve_fit(f, deltas, c11min12_ens)
    C11minC12 = popt[0] / atoms.get_volume()

    # Deriving C11 and C12
    C12 = 1/3 * (3 * B - C11minC12)
    C11 = C12 + C11minC12
    anis = 2 * C44 / (C11 - C12)

    results = {
        'constants': {
            'B': float(B / GPa),
            'C11minC12': float(C11minC12 / GPa),
            'C11': float(C11 / GPa),
            'C12': float(C12 / GPa),
            'C44': float(C44 / GPa),
            'Anisotropy': float(anis),
            'vol': float(atoms.get_volume()),
        },
    }

    return results # Always return a dictionary! Use {} if necessary
