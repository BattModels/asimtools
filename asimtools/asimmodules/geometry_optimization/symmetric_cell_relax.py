#!/usr/bin/env python
'''
Relaxes structure using ASE while retaining symmetry. Requires ASE>=3.23

Author: mkphuthi@github.com
'''
from typing import Dict, Optional
import ase.optimize
from ase.constraints import ExpCellFilter
from ase.spacegroup.symmetrize import FixSymmetry
from ase.io.trajectory import Trajectory
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms

def symmetric_cell_relax(
    calc_id: str,
    image: Dict,
    optimizer: str = 'BFGS',
    fmax: float = 0.003, #Roughly 0.48GPa
    optimizer_args: Optional[Dict] = None,
    fixsymmetry_args: Optional[Dict] = None,
    expcellfilter_args: Optional[Dict] = None,
) -> Dict:
    """Relaxes cell (and atoms) using ase.constraints.ExpCellFilter while retaining symmetry

    :param calc_id: calc_id specification
    :type calc_id: str
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param optimizer: Any optimizer from ase.optimize, defaults to 'BFGS'
    :type optimizer: str, optional
    :param fmax: fmax in optimizer, defaults to 0.003 which is ~0.48GPa
    :type fmax: float, optional
    :param fixsymmetry_args: arguments for ase.constraints.FixSymmetry, defaults to None
    :type fixsymmetry_args: Optional[Dict], optional
    :param expcellfilter_args: arguments for ase.constraints.ExpCellFilter, defaults to None
    :type expcellfilter_args: Optional[Dict], optional
    :return: Results including relaxed structure energy
    :rtype: Dict
    """
    if fixsymmetry_args is None:
        fixsymmetry_args = {}
    if expcellfilter_args is None:
        expcellfilter_args = {}
    if optimizer_args is None:
        optimizer_args = {}

    calc = load_calc(calc_id)
    atoms = get_atoms(**image)
    atoms.set_calculator(calc)

    atoms.set_constraint(FixSymmetry(atoms, **fixsymmetry_args))
    ecf = ExpCellFilter(atoms, **expcellfilter_args)

    dyn = getattr(ase.optimize, optimizer)(ecf)
    traj_file = 'cell_relax.traj'
    traj = Trajectory(
        traj_file,
        'a',
        atoms,
        properties=['energy', 'forces', 'stress'],
    )
    dyn.attach(traj)
    try:
        dyn.run(fmax=fmax)
    except Exception:
        print('Failed to optimize atoms')
        raise

    image_file = 'image_output.xyz'
    atoms.write(image_file, format='extxyz')

    energy = float(atoms.get_potential_energy())
    final_fmax = float(atoms.get_stress().max())

    results = {
        'energy': energy,
        'final_fmax': final_fmax,
        'files':{
            'image': image_file,
            'traj': traj_file,
        }
    }
    return results
