#!/usr/bin/env python
'''
Calculates the monovacancy formation energy from a bulk structure

Author: mkphuthi@github.com
'''

from typing import Dict, Sequence, Optional
import logging
from asimtools.calculators import load_calc
from asimtools.asimmodules.geometry_optimization.atom_relax import atom_relax
from asimtools.asimmodules.geometry_optimization.optimize import optimize
from asimtools.utils import (
    get_atoms,
)

def vacancy_formation_energy(
    calc_id: str,
    image: Dict,
    vacancy_index: int = 0,
    atom_relax_args: Optional[Dict] = None,
    optimize_args: Optional[Dict] = None,
    repeat: Optional[Sequence] = (1,1,1)
) -> Dict:
    """Calculates the monovacancy formation energy from a bulk structure

    :param calc_id: calc_id specification
    :type calc_id: str
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param vacancy_index: Index of atom to remove, defaults to 0
    :type vacancy_index: int, optional
    :param atom_relax_args: Args to pass to
        :func:`asimtools.asimmodules.geometry_optimization.atom_relax.atom_relax`,
        this applies if you do not want to optimize the cell at all
        defaults to None
    :type atom_relax_args: Optional[Dict], optional
    :param optimize_args: Args to pass to
        :func:`asimtools.asimmodules.geometry_optimization.optimize.optimize`,
        this applies if you want to optimize both the cell and positions
        defaults to None
    :type optimize_args: Optional[Dict], optional
    :param repeat: Repeat to apply to image to create supercell, defaults to
        None
    :type repeat: Optional[Sequence], optional
    :return: Empty dictionary
    :rtype: Dict
    """    

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
    
    if atom_relax_args is not None:
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
    else:
        try:
            relax_results = optimize(
                    calc_id=calc_id,
                    image={'atoms': vacant},
                    # optimizer=optimize_args.get('optimizer', 'BFGS'),
                    # fmax=atom_relax_args.get('fmax', 0.003),
                    # prefix=f'vacant_optimized',
                    **optimize_args,
                )
            logging.info('Relaxed vacant structure')
        except ValueError:
            logging.error("Failed to calculate energy of bulk structure, "
                        "see calculator output")
            raise

    vacant_e = relax_results['energy']
    vac_form_e = vacant_e - (len(vacant) / len(bulk) * bulk_e)

    res_dict = {}
    res_dict['bulk_e'] = float(bulk_e)
    res_dict['vacant_e'] = float(vacant_e)
    res_dict['bulk_n'] = len(bulk)
    res_dict['vfe'] = float(vac_form_e)

    results = {'vfe': res_dict}
    return results
