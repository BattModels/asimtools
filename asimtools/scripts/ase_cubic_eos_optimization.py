'''
Calculate the equation of state, fit it and extract the equilibrium 
volume, energy and bulk modulus. This script is based on the example provided 
on the ASE website: https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html

author: mkphuthi@github.com
'''

from typing import Dict, Optional
import logging
import matplotlib.pyplot as plt
import pandas as pd
from ase.eos import calculate_eos
from ase.io import read
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms

def ase_cubic_eos_optimization(
    image: Dict,
    calc_id: str,
    npoints: Optional[int] = 5,
    eos_string: Optional[str] = 'sj',
    eps: Optional[float] = 0.04,
    plot: Optional[bool] = True,
) -> Dict:
    """Calculate the equation of state, fit it and extract the equilibrium 
    volume, energy and bulk modulus

    :param image: image to use
    :type image: Dict
    :param calc_id: calculator id
    :type calc_id: str
    :param db_file: ASE database in which to store results
    :type db_file: bulk.db
    :return: results including equilibrium volume, energy and bulk modulus
    :rtype: Dict
    """
    calc = load_calc(calc_id)
    atoms = get_atoms(**image)
    atoms.calc = calc
    v_init = atoms.get_volume()
    traj_file = 'eos.traj'
    eos = calculate_eos(atoms, trajectory=traj_file, eps=eps, npoints=npoints)
    eos.eos_string = eos_string
    v0, e0, B = eos.fit()  # find minimum
    logging.info('Successfully fit EOS')
    # Do one more calculation at the minimum and write to database:
    new_cell = atoms.cell.copy()
    new_cell *= (v0 / atoms.get_volume())**(1 / 3)
    atoms.set_cell(new_cell, scale_atoms=True)
    atoms.get_potential_energy()
    atoms.info['B'] = B
    atoms.info['v0'] = v0
    atoms.write('eos_image_output.xyz')
    
    if plot:
        eos.plot()
        plt.savefig('eos.png')

    traj = read(traj_file, index=':')

    eos_dict = {
        'energies': [],
        'volumes': [],
        'volume_scale_factors': [],
    }
    for struct in traj:
        eos_dict['energies'].append(struct.get_potential_energy())
        eos_dict['volumes'].append(struct.get_volume())
        eos_dict['volume_scale_factors'].append(struct.get_volume() / v_init)

    eos_df = pd.DataFrame(eos_dict)
    eos_df.to_csv('eos.csv')

    results = {
        'v0': float(v0),
        'e0': float(e0),
        'B': float(B),
        'x0': float(v0 / v_init)
    }
    return results
