#!/usr/bin/env python
'''
Relaxes the given atomic structure using ASE's built-in structure
optimizers
'''

from typing import Dict, Sequence
import logging
import json
import numpy as np
import matplotlib.pyplot as plt
from ase.phonons import Phonons
from asimtools.calculators import load_calc
from asimtools.utils import get_atoms

def ase_phonons(
    calc_id: str,
    image: Dict,
    path: str,
    delta: float = 0.01,
    kpts: Sequence[int] = (20, 20, 20),
    supercell: Sequence[int] = (5,5,5),
) -> Dict:
    """Calculates phonon spectrum and DOS using ASE

    :param calc_id: calc_id specification
    :type calc_id: str
    :param image: Image specification, see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param path: Path in BZ for plot
    :type path: str
    :param delta: delta in ASE phonons, defaults to 0.01
    :type delta: float, optional
    :param kpts: kpts in ASE phonons, defaults to (20, 20, 20)
    :type kpts: Sequence[int], optional
    :param supercell: repeat of image to use as supercell, defaults to (5,5,5)
    :type supercell: Sequence[int], optional
    :return: Empty dictionary
    :rtype: Dict
    """
    
    atoms = get_atoms(**image)
    calc = load_calc(calc_id)

    ph = Phonons(atoms, calc, supercell=supercell, delta=delta)
    try:
        ph.run()
    except ValueError:
        logging.error('Failed to run phonon calculation, check calculator')
        raise

    ph.read(acoustic=True)
    ph.clean()

    npoints = 100
    path = atoms.cell.bandpath(path, npoints=npoints)
    bs = ph.get_band_structure(path)

    bs.write('bs.json')

    dos = ph.get_dos(kpts=kpts).sample_grid(npts=100, width=1e-3)
    # dos_data = np.transpose(
    #     np.vstack([dos.get_weights(), dos.get_energies()])
    # )

    dos_file = 'dos.json'
    with open(dos_file, 'w', encoding='utf-8') as f:
        json.dump({
            'weights': dos.get_weights().tolist(),
            'energies': dos.get_energies().tolist(),
            }, f, indent=2
        )

    # Plotting
    fig = plt.figure(1, figsize=(7, 4))
    ax = fig.add_axes([.12, .07, .67, .85])

    emax = np.max(bs.energies) * 1.01
    bs.plot(ax=ax, emin=0.0, emax=emax)

    dosax = fig.add_axes([.8, .07, .17, .85])
    dosax.fill_between(
        dos.get_weights(),
        dos.get_energies(),
        y2=0,
        color='grey',
        edgecolor='k',
        lw=1
    )

    dosax.set_ylim(0, emax)
    dosax.set_yticks([])
    dosax.set_xticks([])
    dosax.set_xlabel("DOS", fontsize=18)
    fig.savefig('phonons.png')
    plt.close(fig)

    return {}
