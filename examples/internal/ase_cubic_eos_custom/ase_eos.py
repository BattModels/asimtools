'''
Calculate the equation of state, fit it and extract the equilibrium 
volume, energy and bulk modulus. This asimmodule is based on the example provided 
on the ASE website: https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html

author: mkphuthi@github.com
'''

from typing import Dict
import logging
from ase.eos import calculate_eos
from ase.db import connect
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms

def ase_eos(
    image: Dict,
    calculator: dict,
    db_file: 'bulk.db'
) -> Dict:
    """Calculate the equation of state, fit it and extract the equilibrium 
    volume, energy and bulk modulus

    :param image: image to use
    :type image: dict
    :param calculator: calculator spec
    :type calculator: dict
    :param db_file: ASE database in which to store results
    :type db_file: bulk.db
    :return: results including equilibrium volume, energy and bulk modulus
    :rtype: dict
    """
    calc = load_calc(calculator)
    db = connect(db_file)
    atoms = get_atoms(**image)
    atoms.calc = calc
    eos = calculate_eos(atoms)
    v, e, B = eos.fit()  # find minimum
    logging.info('Successfully fit EOS')
    # Do one more calculation at the minimum and write to database:
    atoms.cell *= (v / atoms.get_volume())**(1 / 3)
    atoms.get_potential_energy()
    db.write(atoms, bm=B)

    results = {'v': float(v), 'e': float(e), 'B': float(B)}
    return results
