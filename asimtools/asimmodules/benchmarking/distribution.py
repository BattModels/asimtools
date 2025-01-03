#!/usr/bin/env python
'''
Plots the distribution of various properties in a dataset of structures

Author: mkphuthi@github.com

'''
from typing import Dict, Optional
import numpy as np
import matplotlib.pyplot as plt
from ase.units import kg, m as meters
from asimtools.utils import (
    get_images,
)

cm3 = (meters * 100)**3

def distribution(
    images: Dict,
    unit: str = 'eV',
    bins: int = 50,
    log: bool = True,
    remap_keys: Optional[Dict] = None,
    skip_failed: bool = False,
) -> Dict:
    if remap_keys is None:
        remap_keys = {}
    unit_factors = {'meV': 1000, 'eV': 1, 'kcal/mol': 23.0621}
    unit_factor = unit_factors[unit]

    unit_dict = {
        'energy': f'{unit}/atom',
        'forces': f'{unit}'+r'/$\AA$',
        'stress': f'{unit}'+r'/$\AA^3$',
        'volume': r'$\AA^3$/atom',
        'pressure': f'{unit}'+r'/$\AA^3$',
        'enthalpy': f'{unit}/atom',
        'density': f'g/cm^3',
        'mass': 'amu',
        'natoms': 'atoms',
    }
    images = get_images(**images)
    results = {prop: [] for prop in unit_dict}
    for i, atoms in enumerate(images):
        include = True
        results['natoms'].append(len(atoms))
        if remap_keys.get('energy', False):
            energy = atoms.info[remap_keys['energy']]
        else:
            energy = atoms.get_potential_energy()
        results['energy'].append(energy)
        if remap_keys.get('forces', False):
            forces = atoms.arrays[remap_keys['forces']]
        else:
            forces = atoms.get_forces()
        results['forces'].extend(
            list(np.array(forces).flatten())
        )
        results['volume'].append(atoms.get_volume())
        if remap_keys.get('stress', False):
            stress = atoms.arrays[remap_keys['stress']]
        elif remap_keys.get('virial', False):
            try:
                stress = atoms.info[remap_keys['virial']] / atoms.get_volume()
            except KeyError:
                print('idx:', i, atoms.info, atoms.arrays)
        else:
            stress = atoms.get_stress(voigt=True)
        results['stress'].extend(
            list(np.array(stress)) * unit_factor
        )
        results['pressure'].append(-np.sum(stress[:3])/3)
        mass = np.sum(atoms.get_masses())
        results['mass'].append(mass)
    
    for prop in unit_dict:
        results[prop] = np.array(results[prop])

    results['density'] = (
        (results['mass'] * kg * 1000) / (results['volume'] / cm3)
    )
    results['enthalpy'] = (
        results['energy'] + results['pressure'] * results['volume']
    )
    for prop in ['energy', 'volume', 'enthalpy']:
        results[prop] = results[prop] / results['natoms']
    
    for prop in ['forces', 'stress', 'pressure', 'energy', 'enthalpy']:
        results[prop] = results[prop] * unit_factor

    for prop in unit_dict:
        with open(f'summary.txt', 'a+') as f:
            f.write(f'{prop} distribution\n')
            f.write(f'Num. values: {len(results[prop])}\n')
            f.write(f'Mean: {np.mean(results[prop])} {unit_dict[prop]}\n')
            f.write(f'Std: {np.std(results[prop])} {unit_dict[prop]}\n')
            f.write(f'Min: {np.min(results[prop])} {unit_dict[prop]}\n')
            f.write(f'Max: {np.max(results[prop])} {unit_dict[prop]}\n')
            f.write('+++++++++++++++++++++++++++++++++++++++++++++++\n')
        fig = plt.figure()
        plt.hist(results[prop], bins=bins, density=True)
        plt.xlabel(f'{prop} [{unit_dict[prop]}]')
        plt.ylabel('Frequency')
        if log:
            plt.yscale('log')
        plt.title(f'{prop} distribution')
        plt.savefig(f'{prop}_distribution.png')
        plt.close(fig)
    return {}
