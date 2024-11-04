#!/usr/bin/env python
'''
Plots the distribution of various properties in a dataset of structures

Author: mkphuthi@github.com

'''
from typing import Dict, List, TypeVar, Sequence
import numpy as np
import matplotlib.pyplot as plt
from ase.units import kg, m as meters
from asimtools.utils import (
    get_images,
)

cm3 = (meters * 100)**3

def distribution(
    images: Dict,
    unit: str = 'meV',
    bins: int = 50,
) -> Dict:
    unit_factors = {'meV': 1000, 'eV': 1, 'kcal/mol': 23.0621}
    unit_factor = unit_factors[unit]
    units = {
        'energy': f'{unit}/atom',
        'forces': f'{unit}'+r'/$\AA$',
        'stress': f'{unit}'+r'/$\AA^3$',
    }
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
        results['natoms'].append(len(atoms))
        results['energy'].append(atoms.get_potential_energy())
        results['forces'].extend(
            list(np.array(atoms.get_forces()).flatten())
        )
        results['volume'].append(atoms.get_volume())
        stress = atoms.get_stress(voigt=True)
        results['stress'].extend(
            list(np.array(stress))
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
        print(prop, results[prop])
        results[prop] = results[prop] / results['natoms']



    for prop in unit_dict:
        fig = plt.figure()
        plt.hist(results[prop], bins=bins, density=True)
        plt.xlabel(f'{prop} [{unit_dict[prop]}]')
        plt.ylabel('Frequency')
        plt.title(f'{prop} distribution')
        plt.savefig(f'{prop}_distribution.png')
        plt.close(fig)
    return {}
