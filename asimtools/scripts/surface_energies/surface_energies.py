#!/usr/bin/env python
'''
Describe the script briefly here. If this is a script that runs multiple steps,
describe it here using reStructuredText to generate autodocs

Cite the papers where the method/script was first introduced here as well

Author: mkphuthi@github.com
'''

from typing import Dict, Sequence, Union, Optional
from copy import copy
import logging
import numpy as np
from pymatgen.core.surface import generate_all_slabs
from pymatgen.io.ase import AseAtomsAdaptor as AAA
from asimtools.calculators import load_calc
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
    return surf_en, slab, converged

def surface_energies(
    image: Dict,
    calc_id: str,
    millers: Union[str,Sequence] = 'all',
    generate_all_slabs_args: Optional[Dict] = None,
) -> Dict:
    '''
    Calculates surface energies of slabs defined by args specified for
    pymatgen.core.surface.generate_all_slabs()
    '''

    calc = load_calc(calc_id)
    bulk = get_atoms(**image)
    bulk.set_calculator(calc)

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

    logging.info('Generated %s distinct slabs', len(slabs))
    n_bulk = len(bulk)
    bulk_e_per_atom = bulk.get_potential_energy() / n_bulk
    bulk.write('bulk.xyz')

    slab_dict = {}
    for s, slab in enumerate(slabs):
        big_slab = slab.copy()
        mindex = slab.miller_index
        miller = f'{mindex[0]}{mindex[1]}{mindex[2]}'
        if millers == 'all' or miller in millers:
            logging.info('Calculating for %s', miller)
            atoms = AAA.get_atoms(big_slab)
            atoms.write(f'{miller}.xyz')

            assert np.allclose(atoms.pbc, (True, True, True)), \
                f'Check pbcs for {miller}: {atoms.pbc}'
            assert atoms.cell[2][2] > pymargs['min_vacuum_size'] + \
                pymargs['min_slab_size'], \
                f'Check layer number and vacuum for {miller}'
            assert miller not in slab_dict, \
                f'Multiple terminations for {miller}'

            slab_dict[miller] = {}
            surf_en, slab_atoms, converged = get_surface_energy(
                atoms, load_calc(calc_id), bulk_e_per_atom
            )

            if converged:
                slab_dict[miller]['surf_energy'] = float(surf_en)
                slab_dict[miller]['natoms'] = len(atoms)
                slab_atoms.write(f'{miller}.xyz')

    results = slab_dict

    return results # Always return a dictionary! Use {} if necessary
