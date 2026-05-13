#!/usr/bin/env python
'''
Calculates surface energies of slabs defined by args specified for
pymatgen.core.surface.generate_all_slabs()

Author: mkphuthi@github.com
'''

from typing import Dict, Sequence, Union, Optional
from copy import copy
import logging
import numpy as np
from pymatgen.core.surface import generate_all_slabs
from pymatgen.io.ase import AseAtomsAdaptor as AAA
from asimtools.calculators import load_calc
from asimtools.asimmodules.geometry_optimization.atom_relax import atom_relax
from asimtools.asimmodules.singlepoint import singlepoint
from asimtools.utils import (
    get_atoms,
)

def get_surface_energy(slab, calc, bulk_e_per_atom):
    a1 = slab.get_cell()[0]
    a2 = slab.get_cell()[1]
    area = np.abs(np.cross(a1, a2)[2])
    slab.calc = calc
    converged = True
    try:
        slab_en = slab.get_potential_energy()
    except:
        converged = False
    if converged:
        nslab = len(slab)
        surf_en = (slab_en - nslab * bulk_e_per_atom) / (2 * area)
    else:
        surf_en = None
    return converged, surf_en, slab_en, area

def surface_energies(
    image: Dict,
    calculator: Dict = None,
    millers: Union[str,Sequence] = 'all',
    atom_relax_args: Optional[Dict] = None,
    generate_all_slabs_args: Optional[Dict] = None,
) -> Dict:
    """Calculates surface energies of slabs defined by args specified for
    pymatgen.core.surface.generate_all_slabs()

    :param calculator: Calculator specification, see :func:`asimtools.calculators.load_calc`
    :type calculator: Dict
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param millers: List of miller indices to consider in the form 'xyz',
        defaults to 'all'
    :type millers: Union[str,Sequence], optional
    :param atom_relax_args: Args to pass to
        :func:`asimtools.asimmodules.geometry_optimization.atom_relax.atom_relax`,
        defaults to None
    :type atom_relax_args: Optional[Dict], optional
    :param generate_all_slabs_args: Args to pass to
        :func:`pymatgen.core.surface.generate_all_slabs`, defaults to None
    :type generate_all_slabs_args: Optional[Dict], optional
    :return: _description_
    :rtype: Dict
    """

    calc = load_calc(calculator=calculator)
    bulk = get_atoms(**image)
    bulk.calc = calc

    default_pymatgen_kwargs = {
        'max_index': 3,
        'min_slab_size': 14,
        'min_vacuum_size': 10,
        'max_normal_search': True,
        'symmetrize': True,
    }

    pymargs = copy(default_pymatgen_kwargs)
    pymargs.update(generate_all_slabs_args)

    bulk_struct = AAA.get_structure(bulk)
    slabs = generate_all_slabs(
        bulk_struct,
        **pymargs
    )
    print('Generated %s slabs', len(slabs))
    logging.info('Generated %s distinct slabs', len(slabs))
    n_bulk = len(bulk)
    bulk_e_per_atom = bulk.get_potential_energy() / n_bulk
    bulk.write('bulk.xyz')

    slab_dict = {}
    for _, slab in enumerate(slabs):
        big_slab = slab.copy()
        mindex = slab.miller_index
        miller = f'{mindex[0]}{mindex[1]}{mindex[2]}'
        if millers == 'all' or miller in millers:
            slab_idx = 0
            while f'{miller}_{slab_idx}' in slab_dict:
                slab_idx += 1
            slab_name = f'{miller}_{slab_idx}'

            logging.info('Calculating for %s with index %s', miller, slab_idx)
            atoms = AAA.get_atoms(big_slab)
            atoms.write(f'miller-{slab_name}.xyz')
            slab_dict[slab_name] = {}

            if atom_relax_args is not None:
                relax_results = atom_relax(
                    calculator=calculator,
                    image={'atoms': atoms},
                    optimizer=atom_relax_args.get('optimizer', 'BFGS'),
                    properties=('energy','forces'),
                    fmax=atom_relax_args.get('fmax', 0.01),
                    prefix=f'{miller}_relaxed'
                )

                atoms = get_atoms(
                    image_file = relax_results.get('files', {}).get('image')
                )

                assert np.allclose(atoms.pbc, (True, True, True)), \
                    f'Check pbcs for {miller}: {atoms.pbc}'
                assert atoms.cell[2][2] > pymargs['min_vacuum_size'] + \
                    pymargs['min_slab_size'], \
                    f'Check layer number and vacuum for {miller}'

                converged, surf_en, slab_en, area = get_surface_energy(
                    atoms, load_calc(calculator=calculator), bulk_e_per_atom
                )

                if converged:
                    slab_dict[slab_name]['surf_energy'] = float(surf_en)
                    slab_dict[slab_name]['natoms'] = len(atoms)
                    slab_dict[slab_name]['slab_energy'] = float(slab_en)
                    slab_dict['bulk_energy_per_atom'] = float(bulk_e_per_atom)
                    slab_dict[slab_name]['area'] = float(area)
                    atoms.write(f'{slab_name}.xyz')
                else:
                    logging.warning('Slab %s not converged', miller)

    assert len(slab_dict) > 0, \
        'No slabs generated. Check your args for miller indices'
    results = {
        'surface_energies': slab_dict,
    }

    return results # Always return a dictionary! Use {} if necessary
