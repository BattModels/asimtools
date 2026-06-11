'''
Shared helpers for the geometry_optimization asimmodules.

This module is intentionally not an asimmodule itself (it has no function
matching its file name), it only provides the common optimization backend so
that trajectory handling, convergence checking and results are consistent
across the atom and cell relaxation asimmodules.

Author: mkphuthi@github.com
'''
import logging
from typing import Dict, Optional, Sequence, TypeVar
import numpy as np
import ase.optimize
from ase.io.trajectory import Trajectory
from asimtools.utils import get_logger, write_atoms

Atoms = TypeVar('Atoms')

def run_optimization(
    atoms: Atoms,
    optimizer: str = 'LBFGS',
    fmax: float = 0.02,
    steps: int = 1000,
    optimize_obj: Optional[object] = None,
    prefix: Optional[str] = None,
    label: str = 'optimize',
    properties: Sequence[str] = ('energy', 'forces', 'stress'),
    optimizer_args: Optional[Dict] = None,
    report_stress: bool = False,
    logger: Optional[logging.Logger] = None,
) -> Dict:
    """Run an ASE optimizer with a standard trajectory and results workflow.

    Shared backend for the geometry_optimization asimmodules so that
    trajectory handling, convergence checking, output files and the results
    dictionary are consistent across atom and cell relaxations.

    :param atoms: Atoms object that energy/forces/stress and the trajectory are
        read from. The calculator must already be attached.
    :type atoms: Atoms
    :param optimizer: Any optimizer class name from ase.optimize, defaults to
        'LBFGS'
    :type optimizer: str, optional
    :param fmax: Convergence threshold passed to the optimizer, defaults to 0.02
    :type fmax: float, optional
    :param steps: Maximum number of optimizer steps, defaults to 1000
    :type steps: int, optional
    :param optimize_obj: Object the optimizer acts on, e.g. an ase.filters
        filter wrapping ``atoms``. Defaults to ``atoms`` itself.
    :type optimize_obj: object, optional
    :param prefix: Prefix for output files, defaults to None
    :type prefix: str, optional
    :param label: Base name used for the trajectory file, defaults to 'optimize'
    :type label: str, optional
    :param properties: Properties stored in the trajectory, defaults to
        ('energy', 'forces', 'stress')
    :type properties: Sequence[str], optional
    :param optimizer_args: Extra keyword arguments for the optimizer class,
        defaults to None
    :type optimizer_args: Dict, optional
    :param report_stress: Whether to include the maximum absolute stress
        component ``final_smax`` in the results, defaults to False
    :type report_stress: bool, optional
    :param logger: Logger to report failures to, defaults to a new logger
    :type logger: logging.Logger, optional
    :raises ValueError: If the optimization does not converge below fmax
    :return: Results including energy, final fmax (and smax), convergence flag
        and output files
    :rtype: Dict
    """
    if optimizer_args is None:
        optimizer_args = {}
    if optimize_obj is None:
        optimize_obj = atoms
    if logger is None:
        logger = get_logger()

    prefix = prefix + '_' if prefix else ''

    traj_file = f'{prefix}{label}.traj'
    dyn = getattr(ase.optimize, optimizer)(optimize_obj, **optimizer_args)
    traj = Trajectory(
        traj_file,
        'a',  # Keep frames from multiple relaxations in the same workdir
        atoms,
        properties=list(properties),
    )
    dyn.attach(traj)
    try:
        converged = bool(dyn.run(fmax=fmax, steps=steps))
    except Exception:
        logger.error('Failed to run optimizer %s', optimizer)
        raise
    finally:
        traj.close()

    energy = float(atoms.get_potential_energy())
    final_fmax = float(np.sqrt((atoms.get_forces() ** 2).sum(axis=1).max()))

    results = {
        'energy': energy,
        'final_fmax': final_fmax,
        'converged': converged,
    }
    if report_stress:
        results['final_smax'] = float(np.abs(atoms.get_stress()).max())

    if not converged:
        logger.error('Failed to relax below fmax=%s', fmax)
        raise ValueError(
            f'Optimization did not converge, final fmax={final_fmax}'
        )

    image_file = f'{prefix}image_output.xyz'
    write_atoms(image_file, atoms, fmt='extxyz')
    results['files'] = {'image': image_file, 'traj': traj_file}
    return results
