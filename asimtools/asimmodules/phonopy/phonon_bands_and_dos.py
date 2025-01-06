from typing import Dict, Tuple, Sequence, Optional, Union
from glob import glob
from pathlib import Path
import os
import numpy as np
from numpy.typing import ArrayLike
# import matplotlib.pyplot as plt
# from ase.io import read
# import phonopy
# from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
# from asimtools.calculators import load_calc
# from asimtools.asimmodules.workflows.image_array import image_array
# from asimtools.utils import get_str_btn, get_images
from asimtools.job import UnitJob

def phonon_bands_and_dos(
    image: Dict,
    calc_id: str,
    calc_env_id: str,
    process_env_id: str,
    supercell: ArrayLike = [10,10,10],
    distance: Optional[float] = 0.02,
    phonopy_save_path: Optional[str] = 'phonopy_save.yaml',
    paths: Optional[ArrayLike] = None, 
    labels: Optional[ArrayLike] = None,
    use_seekpath: Optional[bool] = True,
    npoints: Optional[int] = 51,
    mesh: Optional[Union[ArrayLike,float]] = [20, 20, 20],
) -> Dict:
    """Workflow to calculate phonon bands and density of states using phonopy.

    :param image: Image specification. See :ref:`asimtools.utils.get_image`.
    :type image: Dict
    :param calc_id: calc_id of the calculator to use.
    :type calc_id: str
    :param calc_env_id: env_id of the calculator to use.
    :type calc_env_id: str
    :param process_env_id: env_id for pre- and post-processing.
    :type process_env_id: str
    :param supercell: supercell of specified image, defaults to [10,10,10]
    :type supercell: ArrayLike, optional
    :param distance: Distance in phonopy, defaults to 0.02
    :type distance: float, optional
    :param phonopy_save_path: Where to put the phonopy yaml, defaults to 
        'phonopy_save.yaml'
    :type phonopy_save_path: str, optional
    :param paths: Path through BZ for bands plot. Has shape Np*Nik*3 where Np 
        is the number of disconnected paths, Nik is the number of points in
        the i-th path and the last dimension is the coordinates in kspace,
        e.g. [[(0,0,0),(0,0,1)],[(0,0,0),(0.5,0.5,0.0),(0.25,0.25,0.25)]],
        has two disconnected paths with 2 and 3 points respectively,
        defaults to None
    :type paths: Optional[Sequence[Sequence[Sequence]]], optional
    :param labels: Labels for each point provided in paths, defaults to None
    :type labels: Optional[ArrayLike], optional
    :param use_seekpath: Automatically generate a suggested path, defaults to 
        True
    :type use_seekpath: Optional[bool], optional
    :param npoints: Number of points for band path, defaults to 51
    :type npoints: int, optional
    :param mesh: Size of kpoint mesh for DOS either as integer or array of 
        length 3, defaults to [20, 20, 20]
    :type mesh: Union[ArrayLike,float], optional
    :return: result
    :rtype: Dict
    """
    phonopy_save_path = str(Path(phonopy_save_path).resolve())
    sim_input = {
        'asimmodule': 'workflows.chained',
        'args': {
            'steps': {
                'step-0': {
                    'asimmodule': 'phonopy.generate_phonopy_displacements',
                    'env_id': process_env_id,
                    'args': {
                        'image': image,
                        'supercell': supercell,
                        'distance': distance,
                        'phonopy_save_path': phonopy_save_path,
                    },
                },
                'step-1': {
                    'asimmodule': 'phonopy.forces',
                    'env_id': calc_env_id,
                    'args': {
                        'images': {
                            'pattern': '../step-0/supercell-*',
                            'format': 'vasp',
                        },
                        'calc_id': calc_id,
                    },
                },
                'step-2': {
                    'asimmodule': 'phonopy.phonon_bands_and_dos_from_forces',
                    'env_id': process_env_id,
                    'args': {
                        'supercell_image_file_pattern': '../step-1/id-*/image_output.xyz',
                        'phonopy_save_path': str(Path('../step-0') / phonopy_save_path),
                        'paths': paths,
                        'labels': labels,
                        'use_seekpath': use_seekpath,
                        'npoints': npoints,
                        'mesh': mesh,
                    },
                },
            },
        },
    }
    uj = UnitJob(sim_input)
    uj.submit()

    return {}
