from pathlib import Path
from typing import Dict, Optional, Sequence
from asimtools.job import UnitJob

def full_qha(
    image: Dict,
    calc_id: str,
    phonopy_save_path: str,
    calc_env_id: Optional[str] = None,
    process_env_id: Optional[str] = None,
    ase_cubic_eos_args: Optional[Dict] = None,
    supercell: Sequence = [5,5,5],
    t_max: float = 1000,
    pressure: Optional[float] = None,
) -> Dict:
    """Perform a full Quasiharmonic Approximation and predict thermal 
    properties of a given structure. Calculated properties included 
    vibrational free energy heat capacity, thermal expansion etc.

    :param image: Image specification, see :func:`asimtools.utils.get_atoms` 
    :type image: Dict
    :param calc_id: calc_id
    :type calc_id: str
    :param phonopy_save_path: Path where phonopy save yaml is saved, this file
        is important to keep for easier postprocess/analsyis
    :type phonopy_save_path: str
    :param calc_env_id: env_id for running calculator, defaults to None
    :type calc_env_id: Optional[str], optional
    :param process_env_id: env_id for running light pre/post-processing, 
        defaults to None
    :type process_env_id: Optional[str], optional
    :param ase_cubic_eos_args: arguments to pass to 
        :func:`asimtools.asimmodules.geometry_optimization.ase_cubic_eos.ase_cubic_eos` , defaults to None
    :type ase_cubic_eos_args: Optional[Dict], optional
    :param supercell: supercell to use for finite difference method, 
        defaults to [5,5,5]
    :type supercell: Sequence, optional
    :param t_max: Max. temperature for thermal properties, defaults to 1000
    :type t_max: float, optional
    :param pressure: Pressure to optimize to, defaults to None
    :type pressure: Optional[float], optional
    :return: Nothing
    :rtype: Dict
    """

    phonopy_save_path = str(Path(phonopy_save_path).resolve())
    ase_cubic_eos_args['image'] = image
    ase_cubic_eos_args['calc_id'] = calc_id
    scales = ase_cubic_eos_args.get('scales', False)
    if scales:
        npoints = len(scales)
    else:
        npoints = ase_cubic_eos_args['npoints']
    sim_input = {
        'asimmodule': 'workflows.chained',
        'env_id': process_env_id,
        'args': {
            'steps': {
                'step-0': {
                    'asimmodule': 'geometry_optimization.ase_cubic_eos_optimization',
                    'env_id': calc_env_id,
                    'args': ase_cubic_eos_args,
                },
                'step-1': {
                    'asimmodule': 'workflows.sim_array',
                    'env_id': process_env_id,
                    'args': {
                        'template_sim_input': {
                            'asimmodule': 'workflows.chained',
                            'env_id': process_env_id,
                            'args': {
                                'steps': {
                                    'step-0': {
                                        'asimmodule': 'phonopy.generate_phonopy_displacements',
                                        'env_id': process_env_id,
                                        'args': {
                                            'image': {
                                                'image_file': '../../step-0/eos.traj',
                                                'index': {}
                                            },
                                            'supercell': supercell,
                                            'distance': 0.02,
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
                                            'phonopy_save_path': phonopy_save_path,
                                            'use_seekpath': True,
                                        },
                                    },
                                    'step-3': {
                                        'asimmodule': 'phonopy.thermal_properties',
                                        'env_id': process_env_id,
                                        'args': {
                                            'phonopy_save_path': phonopy_save_path,
                                            't_max': t_max,
                                            'suffix': {}
                                        },
                                    },
                                },
                            },
                        },
                        'key_sequence': ['args', 'steps', 'step-3', 'args', 'suffix'],
                        'array_values': [f'{i:02d}' for i in range(npoints)],
                        'labels': 'values',
                        'secondary_key_sequences': [['args', 'steps', 'step-0', 'args', 'image', 'index']],
                        'secondary_array_values': [[i for i in range(npoints)]],
                    },
                },
                'step-2': {
                    'asimmodule': 'phonopy.qha_properties',
                    'env_id': process_env_id,
                    'args': {
                        'ev_csv': '../step-0/eos.csv',
                        'phonopy_save_path': phonopy_save_path,
                        'thermal_properties_file_pattern': '../step-1/id-*/step-3/thermal_properties-*.yaml',
                        't_max': t_max - 10,
                        'pressure': pressure,
                    },
                },
            },
        },
    }

    uj = UnitJob(sim_input)
    uj.submit()

    return {}
