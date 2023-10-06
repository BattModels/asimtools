#!/usr/bin/env python
'''
Describe the script briefly here. If this is a script that runs multiple steps,
describe it here using reStructuredText to generate autodocs

Cite the papers where the method/script was first introduced here as well

Author: mkphuthi@github.com
'''

from typing import Dict, Optional, Sequence
import numpy as np
import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory
from ase.units import GPa, fs
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.md.npt import NPT
from asimtools.calculators import load_calc
from asimtools.utils import (
    get_atoms,
    get_images,
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
        timestep,
        externalstress=0,
        traj_file=None,
        ttime=25*fs,
        pfactor=(75*fs)**2 * 14*GPa,  #ptime=75*fs
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

def plot_thermo(
    images: Dict,
    props: Sequence = ['epot', 'temp', 'ekin', 'etot', 'press'],
    prefix: str = 'ase_md',
):
    atoms = get_images(**images)
    prop_dict = {prop: [] for prop in props}

    for a, atoms in enumerate(images):
        epot = atoms.get_potential_energy()
        ekin = atoms.get_kinetic_energy()
        etot = epot+ekin
        T = atoms.get_temperature()
        prop_dict['epot'].append(epot)
        prop_dict['ekin'].append(ekin)
        prop_dict['etot'].append(etot)
        prop_dict['temp'].append(T)
        if 'press' in props:
            pressure = np.mean(atoms.get_stress(include_ideal_gas=True)[:3])
            prop_dict['press'].append(pressure)

    stride = images.get('index', ':')

    steps = np.arange(len(prop_dict['epots'])) * stride
    for _, prop in enumerate(props):
        _, ax = plt.subplots()
        ax.plot(steps, prop_dict[prop])
        ax.set_xlabel('Step')
        ax.set_ylabel(prop)
        plt.savefig(f'{prefix}_{prop}.png')

    plt.close(fig='all')

def ase_md(
    calc_id: str,
    image: Dict,
    timestep: float,
    temp: float,
    dynamics: str = 'npt',
    friction: float = 1e-2,
    ttime: float = 25*fs,
    pfactor: float = None, #(75*fs)**2 * 14*GPa, #14 is bulk modulus of material i.e. Li
    nsteps: int = 100,
    externalstress: float = 0,
    plot: bool = True,
) -> Dict:
    '''
    Script does xyz specifically
    '''

    calc = load_calc(calc_id)
    atoms = get_atoms(**image)
    atoms.set_calculator(calc)

    if dynamics == 'langevin':
        atoms, _ = langevin_nvt(
            atoms,
            temp,
            nsteps,
            traj_file='output.traj',
            timestep=timestep,
            friction=friction,
        )
    elif dynamics == 'npt':
        atoms, _ = npt(
            atoms,
            temp,
            nsteps,
            traj_file='output.traj',
            timestep=timestep,
            pfactor=pfactor,
            externalstress=externalstress,
            ttime=ttime,
        )

    if plot:
        plot_thermo({'input_file': 'output.traj'})

    results = {}
    return results
