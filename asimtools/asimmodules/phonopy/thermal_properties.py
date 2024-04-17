from typing import Union, Dict, Optional
import numpy as np
from numpy.typing import ArrayLike
import matplotlib.pyplot as plt
import pandas as pd
import yaml
import phonopy

def thermal_properties(
    phonopy_save_path: str,
    mesh: Union[ArrayLike,float] = [20,20,20],
    t_step: float = 10.0,
    t_min: float = 0.0,
    t_max: float = 1000.0,
    run_mesh_kwargs: Optional[Dict] = None,
    run_thermal_properties_kwargs: Optional[Dict] = None,
    suffix: str = '',
) -> Dict:

    if run_mesh_kwargs is None:
        run_mesh_kwargs = {}
    if run_thermal_properties_kwargs is None:
        run_thermal_properties_kwargs = {}

    phonon = phonopy.load(phonopy_save_path)
    phonon.run_mesh(mesh, **run_mesh_kwargs)
    phonon.run_thermal_properties(
        t_step=t_step,
        t_max=t_max,
        t_min=t_min,
        **run_thermal_properties_kwargs,
    )
    phonon.save(phonopy_save_path)

    # Write and plot thermal properties
    if suffix != '':
        suffix = '-' + suffix
    tp_label = 'thermal_properties' + suffix
    tp_dict = phonon.get_thermal_properties_dict()
    temperatures = tp_dict['temperatures']
    free_energy = tp_dict['free_energy']
    entropy = tp_dict['entropy']
    heat_capacity = tp_dict['heat_capacity']

    tp_df = pd.DataFrame({
        'Temperature (K)': temperatures,
        'Free Energy (kJ/mol)': free_energy,
        'Entropy (J/K/mol)': entropy,
        'Heat Capacity (J/K/mol)': heat_capacity,
    })
    tp_df.to_csv(f'{tp_label}_table.csv', index=False)
    phonon.plot_thermal_properties()
    plt.savefig(f'{tp_label}.png')

    # Write yaml file needed for QHA
    phonon.write_yaml_thermal_properties(
        f'{tp_label}.yaml'
    )

    return {}
