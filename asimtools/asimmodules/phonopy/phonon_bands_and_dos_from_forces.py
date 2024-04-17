from typing import Dict, Tuple, Sequence, Optional, Union
from glob import glob
import os
import numpy as np
from numpy.typing import ArrayLike
import matplotlib.pyplot as plt
from ase.io import read
import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from asimtools.calculators import load_calc
from asimtools.asimmodules.workflows.image_array import image_array
from asimtools.utils import get_str_btn, get_images

def phonon_bands_and_dos_from_forces(
    supercell_image_file_pattern: str,
    phonopy_save_path: str,
    paths: Optional[Sequence] = None,
    labels: Optional[Sequence] = None,
    use_seekpath: Optional[bool] = False,
    npoints: int = 51,
    mesh: Union[ArrayLike,float] = [20, 20, 20],
) -> Dict:
    assert use_seekpath or paths is not None, \
        'Provide explicit band paths or set use_seekpath=True to automate'
    images = get_images(pattern=supercell_image_file_pattern)
    Ncells = len(images)
    Natoms = len(images[0])
    force_set = np.zeros([Ncells, Natoms, 3])
    for i, atoms in enumerate(images):
        assert len(atoms) == Natoms, 'Images must have the same num. of atoms'
        try:
            forces = atoms.get_forces()
        except Exception('Failed to get forces from image %s', i):
            raise
        force_set[i] = np.array(forces)

    phonon = phonopy.load(phonopy_save_path)
    phonon.forces = force_set
    phonon.produce_force_constants()
    phonon.save(phonopy_save_path, settings={'force_constants': True})

    if use_seekpath:
        phonon.auto_band_structure(write_yaml=True)
    else:
        qpoints, connections = get_band_qpoints_and_path_connections(
            paths, npoints=npoints
        )
        phonon.run_band_structure(
            qpoints, path_connections=connections, labels=labels
        )
        phonon.write_yaml_band_structure()

    phonon.plot_band_structure()
    plt.savefig('phonon_bands.png')

    phonon.auto_total_dos(
        write_dat=True,
        xlabel='Frequency (THz)',
        mesh=mesh
    )
    phonon.plot_total_dos()
    plt.savefig('phonon_dos.png')

    phonon.save(phonopy_save_path)
    return {}
