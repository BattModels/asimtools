#!/usr/bin/env python
'''
Describe the asimmodule briefly here. If this is a asimmodule that runs multiple steps,
describe it here using reStructuredText to generate autodocs

Cite the papers where the method/asimmodule was first introduced here as well

Author: mkphuthi@github.com
'''

from typing import Dict, Sequence, Optional
import logging
from asimtools.calculators import load_calc
from asimtools.asimmodules.geometry_optimization.atom_relax import atom_relax
from asimtools.utils import (
    get_atoms,
)

def vacancy_formation_energy(
    calc_id: str,
    image: Dict,
    vacancy_index: int = 0,
    atom_relax_args: Optional[Dict] = None,
    repeat: Optional[Sequence] = None
) -> Dict:

    calc = load_calc(calc_id)
    bulk = get_atoms(**image).repeat(repeat)
    bulk.set_calculator(calc)

    vacant = bulk.copy()
    del vacant[vacancy_index]

    try:
        bulk_e = calc.get_potential_energy(bulk)
        logging.info('Finished bulk calculation')
    except ValueError:
        logging.error("Failed to calculate energy of bulk structure, "
                      "see calculator output")
        raise

    bulk.write('bulk_image_output.xyz')
    
    try:
        relax_results = atom_relax(
                calc_id=calc_id,
                image={'atoms': vacant},
                optimizer=atom_relax_args.get('optimizer', 'BFGS'),
                properties=('energy','forces'),
                fmax=atom_relax_args.get('fmax', 0.01),
                prefix=f'vacant_relaxed',
            )
        logging.info('Relaxed vacant structure')
    except ValueError:
        logging.error("Failed to calculate energy of bulk structure, "
                      "see calculator output")
        raise

    vacant_e = relax_results['energy']
    vac_form_e = vacant_e - (len(vacant) / len(bulk) * bulk_e)

    res_dict = {}
    res_dict['bulk_e'] = bulk_e
    res_dict['vacant_e'] = vacant_e
    res_dict['bulk_n'] = len(bulk)
    res_dict['vfe'] = vac_form_e

    results = {'vfe': res_dict}
    return results
