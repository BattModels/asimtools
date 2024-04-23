'''
Calculate the equation of state, fit it and extract the equilibrium 
volume, energy and bulk modulus. This script is based on the example provided 
on the ASE website: https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html

author: mkphuthi@github.com
'''

from typing import Dict, Optional, Sequence
import logging
import matplotlib.pyplot as plt
import pandas as pd
from ase.eos import EquationOfState, calculate_eos
from ase.io import read
from ase.io.trajectory import Trajectory
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms

def ase_cubic_eos_optimization(
    calc_id: str,
    image: Dict,
    npoints: Optional[int] = 5,
    eos_string: Optional[str] = 'sj',
    eps: Optional[float] = 0.04,
    scales: Optional[Sequence] = None,
    plot: Optional[bool] = True,
) -> Dict:
    """Generate the energy-volume equation of state (energy calculations not parallelized)

    :param calc_id: calc_id specification
    :type calc_id: str
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param npoints: Number of energy points to calculate, must be >5, defaults to 5
    :type npoints: Optional[int], optional
    :param eos_string: eos_string as specified for :func:`ase.eos.calculate_eos`, defaults to 'sj'
    :type eos_string: Optional[str], optional
    :param eps: eps as sepecified for ase, defaults to 0.04
    :type eps: Optional[float], optional
    :param scales: array of values to scale the unit cell by (not volume)
    :type eps: Optional[float], optional
    :param plot: Whether to plot eos or not, defaults to True
    :type plot: Optional[bool], optional
    :return: Equilibrium energy, volume, bulk modulus and factor by which to scale lattice parameter to get equilibrium  structure
    :rtype: Dict
    """
    calc = load_calc(calc_id)
    atoms = get_atoms(**image)
    atoms.calc = calc
    v_init = atoms.get_volume()
    traj_file = 'eos.traj'
    if scales is not None:
        traj = Trajectory('eos.traj', 'w')
        volumes = []
        energies = []
        images = []
        for scale in scales:
            scaled_atoms = atoms.copy()
            cell = scaled_atoms.cell
            scaled_atoms.set_cell(cell * scale, scale_atoms=True)
            scaled_atoms.set_calculator(calc)
            volumes.append(scaled_atoms.get_volume())
            energies.append(scaled_atoms.get_potential_energy())
            traj.write(scaled_atoms)

        eos = EquationOfState(volumes, energies)
    else:
        eos = calculate_eos(atoms, trajectory=traj_file, eps=eps, npoints=npoints)


    traj = read(traj_file, index=':')
    eos_dict = {
        'energies': [],
        'volumes': [],
        'volume_scale_factors': [],
    }
    for struct in traj:
        eos_dict['energies'].append(struct.get_potential_energy())
        eos_dict['volumes'].append(struct.get_volume())
        v_factor = struct.get_volume() / v_init
        eos_dict['volume_scale_factors'].append(v_factor)

    eos_df = pd.DataFrame(eos_dict)
    eos_df.to_csv('eos.csv')

    eos.eos_string = eos_string
    v0, e0, B = eos.fit()  # find minimum

    logging.info('Successfully fit EOS')
    # Do one more calculation at the minimum and write to disk
    new_cell = atoms.cell.copy()
    x = (v0 / atoms.get_volume())**(1 / 3) # scale cell by this much
    new_cell *= x
    atoms.set_cell(new_cell, scale_atoms=True)
    atoms.get_potential_energy()
    atoms.info['B'] = B
    atoms.info['v0'] = v0
    atoms.info['scale'] = x
    atoms.write('eos_image_output.xyz')

    if plot:
        eos.plot()
        plt.savefig('eos.png')

    results = {
        'v0': float(v0),
        'e0': float(e0),
        'B': float(B),
        'x0': float(v0 / v_init)
    }
    return results
