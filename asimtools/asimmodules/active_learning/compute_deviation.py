from typing import Dict, List, Sequence, Optional
import os
from glob import glob
from natsort import natsorted
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ase.io import write
from asimtools.calculators import load_calc
from asimtools.utils import (
    change_dict_value,
    get_images,
)

def compute_deviation(
    images: Dict,
    template_calc_params: Optional[Dict] = None,
    model_weights_key_sequence: Optional[Sequence] = None,
    model_weights_pattern: Optional[os.PathLike] = None,
    calc_ids: Optional[Sequence] = None,
) -> Dict:
    """Computes variance of properties from a trajectory file

    :param images: Images specification, see :func:`asimtools.utils.get_images`
    :type images: Dict
    :param template_calc_params: Template calc_params, defaults to None
    :type template_calc_params: Optional[Dict]
    :param model_weights_key_sequence: Sequence of keys to change in the
        template calc_params
    :type model_weights_key_sequence: Optional[Sequence]
    :param model_weights_pattern: Pattern of model weights files, defaults to
        None
    :type model_weights_pattern: Optional[os.PathLike]
    :param calc_ids: List of calc_ids to use, if provided, all other arguments
        are ignored, defaults to None
    :type calc_ids: Optional[Sequence]

    """
    properties = ['energy', 'forces', 'stress', 'energy_per_atom']
    if calc_ids is None:
        model_weights_files = natsorted(glob(model_weights_pattern))

        calc_dict = {}
        for i, model_weights_file in enumerate(model_weights_files):
            new_calc_params = change_dict_value(
                template_calc_params,
                model_weights_file,
                key_sequence=model_weights_key_sequence,
                return_copy=True
            )

            calc_dict[f'calc-{i}'] = new_calc_params
    else:
        calc_dict = {calc_id: calc_id for calc_id in calc_ids}

    variances = {prop: {} for prop in properties}

    images = get_images(**images)
    assert len(images) > 0, 'No images found'

    write('images_input.xyz', images)

    prop_dict = {
        prop: {calc_id: [] for calc_id in calc_dict} for prop in properties
    }
    for prop in properties:
        if prop == 'forces':
            prop_dict[prop]['mean_std'] = []
            prop_dict[prop]['max_std'] = []
        else:
            prop_dict[prop]['std'] = []

    for i, atoms in enumerate(images):
        atom_results = {prop: [] for prop in properties}
        for calc_id in calc_dict:
            # Some calculators behave badly if not reloaded unfortunately
            if isinstance(calc_dict[calc_id], str):
                calc = load_calc(calc_id=calc_dict[calc_id])
            else:
                calc = load_calc(calc_params=calc_dict[calc_id].copy())

            atoms.set_calculator(calc)
            energy = atoms.get_potential_energy(atoms)
            energy_per_atom = energy / len(atoms)
            forces = np.linalg.norm(atoms.get_forces(), axis=1)
            stress = -np.sum(
                atoms.get_stress(voigt=True, include_ideal_gas=False)[:3]
            ) / 3

            prop_dict['energy_per_atom'][calc_id].append(energy_per_atom)
            prop_dict['energy'][calc_id].append(energy)
            prop_dict['forces'][calc_id].append(forces)
            prop_dict['stress'][calc_id].append(stress)
            atom_results['energy_per_atom'].append(energy_per_atom)
            atom_results['energy'].append(energy)
            atom_results['forces'].append(forces)
            atom_results['stress'].append(stress)

        prop_dict['energy']['std'].append(np.std(atom_results['energy']))
        prop_dict['energy_per_atom']['std'].append(
            np.std(atom_results['energy_per_atom'])
        )
        prop_dict['forces']['mean_std'].append(
            np.mean(np.std(atom_results['forces'], axis=1))
        )
        prop_dict['forces']['max_std'].append(
            np.max(np.std(atom_results['forces'], axis=1))
        )
        prop_dict['stress']['std'].append(np.std(atom_results['stress']))

    df = pd.DataFrame({
        'energy_std': prop_dict['energy']['std'],
        'energy_per_atom_std': prop_dict['energy_per_atom']['std'],
        'force_mean_std': prop_dict['forces']['mean_std'],
        'force_max_std': prop_dict['forces']['max_std'],
        'stress_std': prop_dict['stress']['std']

    })
    df.to_csv('deviations.csv')

    unit_dict = {
        'energy': 'eV',
        'energy_per_atom': 'eV/atom',
        'forces': 'eV/$\AA$',
        'stress': 'eV/$\AA^3$',
    }
    for prop in properties:
        if prop not in ['forces']:
            # df = pd.DataFrame(prop_dict[prop])
            fig, ax = plt.subplots()
            for calc_id in calc_dict:
                ax.plot(prop_dict[prop][calc_id], label=calc_id)
            ax.set_ylabel(f'{prop} {unit_dict[prop]}')
        else:
            fig, ax = plt.subplots()
            for calc_id in calc_dict:
                ax.plot(
                    np.max(prop_dict[prop][calc_id],axis=1),
                    label=calc_id,
                )
            ax.set_ylabel(f'{prop} max {unit_dict[prop]}')
        ax.set_xlabel('Image index')    
        ax.legend()
        plt.savefig(f'{prop}.png')
        plt.close()

        fig, ax = plt.subplots()
        if prop == 'forces':
            ax.plot(prop_dict[prop]['mean_std'], label='mean_std')
            ax.plot(prop_dict[prop]['max_std'], label='max_std')
            ax.legend()
        else:
            ax.plot(prop_dict[prop]['std'])
        ax.set_xlabel('Image index')
        ax.set_ylabel(f'{prop} std {unit_dict[prop]}')
        plt.savefig(f'{prop}_std.png')
        plt.close()
    
    return {}
    