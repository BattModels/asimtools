#!/usr/bin/env python
'''
Generates a parity plot and collects evaluation statistics comparing energy
and/or forces and/or stress to existing values in the provided dataset. This
asimmodule can work in parallel based on the number of cores specified.

Author: mkphuthi@github.com

'''
from typing import Dict, List, TypeVar, Sequence
from functools import partial
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
from asimtools.calculators import load_calc
from asimtools.utils import (
    write_csv_from_dict,
    get_images,
    parse_slice,
    get_axis_lims,
)

Calculator = TypeVar('Calculator')
def calc_parity_data(
    subset: List,
    calc_id: str,
    properties: Sequence = ('energy', 'forces', 'stress'),
    force_prob: float = 1.0,
) -> Dict:
    """Calculates parity data for each atoms instance in subset

    :param subset: List of atoms instances
    :type subset: List
    :param calc_id: calc_id specification
    :type calc_id: str
    :param properties: Properties to evaluate, choose from "energy", \
        "forces" and "stress", defaults to ('energy', 'forces', 'stress')
    :type properties: List, optional
    :param force_prob: Fraction of forces to consider, may speed up sampling \
        large structures by subsampling forces randomly, defaults to 1.0
    :type force_prob: float, optional
    :return: Dictionary with reference and predicted values for each \
        property
    :rtype: Dict
    """

    res = {prop: [] for prop in properties}
    ervals = []
    epvals = []
    frvals = []
    fpvals = []
    srvals = []
    spvals = []
    for i, atoms in enumerate(subset):
        calc = load_calc(calc_id)
        n_atoms = len(atoms)
        if 'energy' in properties:
            prop = 'energy'
            ervals = np.hstack(
                [ervals, atoms.get_potential_energy()/n_atoms]
            )
            epvals = np.hstack(
                [epvals, float(calc.get_potential_energy(atoms)/n_atoms)]
            )

        if 'forces' in properties:
            if np.random.rand() < force_prob:
                prop = 'forces'
                frvals = np.hstack(
                    [frvals, np.array(atoms.get_forces()).flatten()]
                )

                fpvals = np.hstack(
                    [fpvals, np.array(calc.get_forces(atoms)).flatten()]
                )

        if 'stress' in properties:
            prop = 'stress'
            srvals = np.hstack(
                [srvals, np.array(atoms.get_stress()).flatten()]
            )

            spvals = np.hstack(
                [spvals, np.array(calc.get_stress(atoms)).flatten()]
            )
            res[prop] = {'ref': srvals, 'pred': spvals}

        if i % 20 == 0:
            print(f'Progress {i}/{len(subset)}')

    res['energy'] = {'ref': ervals, 'pred': epvals}
    res['forces'] = {'ref': frvals, 'pred': fpvals}
    res['stress'] = {'ref': srvals, 'pred': spvals}
    return res

def _split_data(data: Sequence, chunks: int) -> List:
    """Splits data into equal chunks and tosses the remainder

    :param data: Sequence of datapoints
    :type data: Sequence
    :param chunks: Number of chunks to split the data into
    :type nprocs: int
    :return: List of split datasets
    :rtype: List
    """
    subsets = []
    subset_len = int(np.ceil(len(data) // chunks))
    for subset_id in range(chunks):
        start = subset_id * subset_len
        stop = (subset_id + 1) * subset_len
        subsets.append(data[start:stop])

    return subsets

def rmse(yhat: Sequence, y: Sequence) -> float:
    """Calculate Root Mean Square Error between to sequences

    :param yhat: Predicted values
    :type yhat: Sequence
    :param y: True values
    :type y: Sequence
    :return: RMSE
    :rtype: float
    """
    return np.sqrt(np.sum(np.square(yhat - y)) / len(y))

def parity(
    images: Dict,
    calc_id: str,
    force_prob: float = 1.0,
    nprocs: int = 1,
    unit: str = 'meV',
    index: str = ':',
    properties: Sequence = ('energy', 'forces', 'stress'),
) -> Dict:
    """Generates a parity plot and collects evaluation statistics comparing energy
    and/or forces and/or stress to existing values in the provided dataset

    :param images: Image specification, see :func:`asimtools.utils.get_images`
    :type images: Dict
    :param calc_id: ID of calculator provided in calc_input or global file
    :type calc_id: str
    :param force_prob: Fraction of forces to consider in force parity, \
        can be used for speeding up large structures by only subsampling\
        randomly, defaults to 1.0
    :type force_prob: float, optional
    :param nprocs: Number of process to parallelize over, must be less than \
        number of datapoints, defaults to 1
    :type nprocs: int, optional
    :param unit: Unit to plot, choose from meV, eV or kcal/mol, \
        defaults to 'meV'
    :type unit: str, optional
    :param index: index of data to use from full dataset, defaults to ':'
    :type index: int, optional
    :param properties: Properties to evaluate, choose from "energy",\
         "forces" and "stress", defaults to ('energy', 'forces', 'stress')
    :type properties: Sequence, optional
    :return: Dictionary of results
    :rtype: Dict
    """

    index = parse_slice(index)
    data = get_images(**images)[index]
    assert len(data) > nprocs, \
        f'Too little samples ({len(data)})  for {nprocs}'

    unit_factors = {'meV': 1000, 'eV': 1, 'kcal/mol': 23.0621}
    unit_factor = unit_factors[unit]
    units = {
        'energy': f'{unit}/atom',
        'forces': f'{unit}'+r'/$\AA$',
        'stress': f'{unit}'+r'/$\AA^3$',
    }

    subsets = _split_data(data, nprocs)
    reses = []
    with Pool(nprocs) as pool:
        reses = pool.map(partial(
            calc_parity_data,
            calc_id=calc_id,
            properties=properties,
            force_prob=force_prob,
            ),
        subsets,
        )

    res = {prop: {'ref': [], 'pred': []} for prop in properties}
    results = {}
    for prop in properties:
        for source in ('ref', 'pred'):
            res[prop][source] = np.hstack(
                [reses[jj][prop][source] for jj in range(len(reses))]
            ) * unit_factor

        write_csv_from_dict(
            f'{prop}_parity.csv',
            res[prop],
            header=f'Units: {unit}',
        )

        rmse_val = rmse(res[prop]['ref'], res[prop]['pred'])
        results[f'{prop}_rmse'] = float(rmse_val)

        _, ax = plt.subplots()
        lims = np.array(get_axis_lims(res[prop]['ref'], res[prop]['pred']))

        ax.plot(lims, lims, 'k-')
        ax.plot(
            res[prop]['ref'],
            res[prop]['pred'],
            'ro',
            label=f'{prop.capitalize()} RMSE={rmse_val:.3g}{units[prop]}',
        )
        ax.set_xticks(np.round(np.linspace(*lims, 5), 1))
        ax.set_yticks(np.round(np.linspace(*lims, 5), 1))
        ax.set_aspect('equal', 'box')
        ax.set_xlabel(f'Reference {prop.capitalize()} ({units[prop]})')
        ax.set_ylabel(f'Predicted {prop.capitalize()} ({units[prop]})')
        ax.legend()
        plt.savefig(f'{prop}_parity.png')

    results['Energy unit'] = unit
    return results
