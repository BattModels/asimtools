#!/usr/bin/env python
'''
Describe the script briefly here. If this is a script that runs multiple steps,
describe it here using reStructuredText to generate autodocs

Cite the papers where the method/script was first introduced here as well

Author: mkphuthi@github.com
'''

from __future__ import annotations
from typing import Dict, Optional, Sequence
import numpy as np
import matplotlib.pyplot as plt
import ase
from ase.io.trajectory import Trajectory
from ase.units import GPa, fs
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.md.npt import NPT
from asimtools.calculators import load_calc
from asimtools.utils import (
    get_atoms,
    get_images,
    expand_wildcards,
)

def langevin_nvt(
    atoms: ase.Atoms,
    temp: float,
    nsteps: int,
    traj_file: str = None,
    friction: float = 1e-2,
    timestep: float = 1*fs,
):
    """Does Langevin dynamics

    :param atoms: atoms object
    :type atoms: ase.Atoms
    :param temp: Temperature in Kelvin
    :type temp: float
    :param nsteps: Number of steps to run
    :type nsteps: int
    :param traj_file: trajectory file name, defaults to None
    :type traj_file: str, optional
    :param friction: friction parameter, defaults to 1e-2
    :type friction: float, optional
    :param timestep: Timestep in ASE time units, defaults to 1*fs
    :type timestep: float, optional
    :return: final atoms object and trajectory object
    :rtype: ase.Atoms, Trajectory
    """
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
    atoms: ase.Atoms,
    temp: float,
    nsteps: int,
    timestep: float,
    externalstress: float = 0,
    traj_file: str = None,
    ttime: float = 25*fs,
    pfactor: Optional[float] = None, # (75*fs)**2 * 14*GPa, #Replace 14 with bulk modulus of material
):
    """Does NPT dynamics

    :param atoms: input atoms object
    :type atoms: ase.Atoms
    :param temp: Temperature in Kelvin
    :type temp: float
    :param nsteps: Number of steps to run
    :type nsteps: int
    :param timestep: Timestep in ASE time units
    :type timestep: float
    :param externalstress: External pressure to apply, defaults to 0
    :type externalstress: float, optional
    :param traj_file: trajectory file name, defaults to None
    :type traj_file: str, optional
    :param ttime: thermostat time constant, defaults to 25*fs
    :type ttime: float, optional
    :return: final atoms object and trajectory object
    :rtype: ase.Atoms, Trajectory
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
    """Plots thermodynamic properties of a trajectory

    :param images: Dictionary of images to read from
    :type images: Dict
    :param props: Properties to plot, defaults to ['epot', 'temp', 'ekin', 'etot', 'press']
    :type props: Sequence, optional
    :param prefix: Prefix for output file names, defaults to 'ase_md'
    :type prefix: str, optional
    """
    images_spec = images.copy()
    images = get_images(**images)
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

    steps = np.arange(len(prop_dict['epot']))
    for _, prop in enumerate(props):
        _, ax = plt.subplots()
        ax.plot(steps, prop_dict[prop])
        ax.set_xlabel('Step')
        ax.set_ylabel(prop)
        plt.savefig(f'{prefix}_{prop}.png')

    plt.close(fig='all')

def ase_md(
    calc_spec: Dict,
    image: Dict,
    timestep: float,
    nsteps: int = 100,
    dynamics: str = 'npt',
    temp: Optional[float] = None,
    friction: Optional[float] = 1e-2,
    ttime: Optional[float] = 25*fs,
    pfactor: Optional[float] = None,
    externalstress: Optional[float] = 0,
    plot: Optional[bool] = True,
) -> Dict:
    """Runs ASE MD simulations. This is only recommended for small systems and
    for testing. For larger systems, use LAMMPS or more purpose-built code

    :param calc_spec: calc specification, any entries specified with asterisk 
        wildcards will be expanded
    :type calc_spec: Dict
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param timestep: Timestep in ASE time units
    :type timestep: float
    :param temp: Temperature in Kelvin
    :type temp: float
    :param dynamics: Type of dynamics to run from 'nvt', 'langevin' and 'nvt',
        defaults to 'npt'
    :type dynamics: str, optional
    :param friction: Friction parameter, defaults to 1e-2
    :type friction: float, optional
    :param ttime: Thermostat time constant, defaults to 25*fs
    :type ttime: float, optional
    :param pfactor: Pressure factor, defaults to None
    :type pfactor: float, optional
    :param nsteps: Number of steps to run, defaults to 100
    :type nsteps: int, optional
    :param externalstress: External stress to apply, defaults to 0
    :type externalstress: float, optional
    :param plot: Whether to plot thermodynamic properties, defaults to True
    :type plot: bool, optional
    :return: Results dictionary
    :rtype: Dict
    """

    calc_spec = expand_wildcards(calc_spec)
    calc = load_calc(**calc_spec)
    atoms = get_atoms(**image)
    atoms.set_calculator(calc)

    assert dynamics in ['nvt', 'langevin', 'npt'], 'Invalid dynamics'
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
        assert pfactor is not None, 'Pressure factor must be provided'
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
    elif dynamics == 'nvt':
        atoms, _ = npt(
            atoms,
            temp,
            nsteps,
            traj_file='output.traj',
            timestep=timestep,
            pfactor=None,
            ttime=ttime,
        )

    if plot:
        plot_thermo(images={'image_file': 'output.traj'})

    results = {}
    return results
