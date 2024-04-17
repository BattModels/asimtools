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
    distance: float = 0.02,
    phonopy_save_path: str = 'phonopy_save.yaml',
    paths: Optional[Sequence] = None,
    labels: Optional[Sequence] = None,
    use_seekpath: Optional[bool] = True,
    npoints: int = 51,
    mesh: Union[ArrayLike,float] = [20, 20, 20],
    t_max: float = 1000,
) -> Dict:

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
