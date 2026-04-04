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
    phonopy_save_path: os.PathLike,
    paths: Optional[Sequence[Sequence[Sequence]]] = None,
    labels: Optional[Sequence] = None,
    use_seekpath: Optional[bool] = False,
    npoints: int = 51,
    mesh: Union[Sequence,float] = [20, 20, 20],
) -> Dict:
    """Calculate phonon bands and DOS from already calculated forces. 
    The images with calculatd forces should correspond to the displacements 
    in the phonopy save file

    :param supercell_image_file_pattern: images where forces are calculated
    :type supercell_image_file_pattern: str
    :param phonopy_save_path: phonopy save file path with corresponding displacements
    :type phonopy_save_path: str
    :param paths: Path through BZ for bands plot, given as a list of 
        disconnected paths of arbitrary length where each path is specified by
        a series of coordinates 
        e.g. [[(0,0,0),(0,0,1)],[(0,0,0),(0.5,0.5,0.0),(0.25,0.25,0.25)]],
        defaults to None
    :type paths: Optional[Sequence[Sequence[Sequence]]], optional
    :param labels: What to label the special points specified in the path,
        defaults to None
    :type labels: Optional[Sequence], optional
    :param use_seekpath: Auto-generate a bandpath using seekpath,
        defaults to False
    :type use_seekpath: Optional[bool], optional
    :param npoints: How many points to sample for each segment, defaults to 51
    :type npoints: int, optional
    :param mesh: Mesh to use for binning DOS, defaults to [20, 20, 20]
    :type mesh: Union[ArrayLike,float], optional
    :return: Results
    :rtype: Dict
    """
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
