#!/usr/bin/env python
'''
Describe the script briefly here. If this is a script that runs multiple steps,
describe it here using reStructuredText to generate autodocs

Cite the papers where the method/script was first introduced here as well

Author: mkphuthi@github.com
'''

from typing import Dict, Optional
import numpy as np
from ase.io.trajectory import Trajectory
from ase.units import GPa, fs
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.md.npt import NPT
from asimtools.calculators import load_calc
from asimtools.utils import (
    get_atoms,
    # get_images,
)

def langevin_nvt(
    atoms, temp, nsteps, traj_file=None, friction=1e-2, timestep=1*fs):
    MaxwellBoltzmannDistribution(atoms, temperature_K=temp)
    dyn = Langevin(
        atoms,
        timestep=timestep,
        temperature_K=temp,
        friction=friction,
    )

    if traj_file is not None:
        traj = Trajectory(
            traj_file,
            atoms=atoms,
            mode='w',
            properties=['energy', 'forces', 'stress']
        )
        dyn.attach(traj.write)
    dyn.run(nsteps)
    return atoms, traj

def npt(
        atoms,
        temp,
        nsteps,
        externalstress=0,
        traj_file=None,
        ttime=25*fs,
        pfactor=(75*fs)**2 * 14*GPa,  #ptime=75*fs
        timestep=1*fs,
    ):
    """Does Nose-Hoover dynamics

    To find timestep from Nose mass (Q), tau^2 = Q / (Ndof * kB * T0)
    :param atoms: _description_
    :type atoms: _type_
    :param temp: _description_
    :type temp: _type_
    :param nsteps: _description_
    :type nsteps: _type_
    :param press: _description_, defaults to 0
    :type press: int, optional
    :param traj_file: _description_, defaults to None
    :type traj_file: _type_, optional
    :param ttime: _description_, defaults to 25*fs
    :type ttime: _type_, optional
    :param pfactor: _description_, defaults to (75*fs)**2*14*GPa
    :type pfactor: _type_, optional
    :param externalstress: _description_, defaults to None
    :type externalstress: _type_, optional
    :return: _description_
    :rtype: _type_
    """
    MaxwellBoltzmannDistribution(atoms, temperature_K=temp)
    dyn = NPT(
        atoms,
        timestep=timestep,
        temperature_K=temp,
        externalstress=externalstress,
        ttime=ttime,
        pfactor=pfactor,
        # mask=np.diag([1, 1, 1]),
    )
    if traj_file is not None:
        traj = Trajectory(
            traj_file,
            atoms=atoms,
            mode='w',
            properties=['energy', 'forces', 'stress'],
        )
        dyn.attach(traj.write)
    dyn.run(nsteps)
    traj = Trajectory(traj_file, 'r')
    return atoms, traj

def ase_md(
    calc_id: str,
    image: Dict,
    dynamics: str = 'npt',
    temp: float = 300.0,
    timestep: float = 1*fs,
    friction: float = 1e-2,
    ttime: float = 25*fs,
    pfactor: float = None, #(75*fs)**2 * 14*GPa, #14 is bulk modulus of material i.e. Li
    nsteps: int = 100,
    externalstress: float = 0,
) -> Dict:
    '''
    Script does xyz specifically
    '''

    calc = load_calc(calc_id)
    atoms = get_atoms(**image)
    atoms.set_calculator(calc)

    if dynamics == 'langevin':
        atoms, traj = langevin_nvt(
            atoms,
            temp,
            nsteps,
            traj_file='output.traj',
            timestep=timestep,
            friction=friction,
        )
    elif dynamics == 'npt':
        atoms, traj = npt(
            atoms,
            temp,
            nsteps,
            traj_file='output.traj',
            timestep=timestep,
            pfactor=pfactor,
            externalstress=externalstress,
            ttime=ttime,
        )

    results = {}
    return results
