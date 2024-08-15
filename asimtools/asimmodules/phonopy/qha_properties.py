"""Example of QHA calculation by Al."""
from typing import Union, Dict, Optional
from glob import glob
import os
import numpy as np
from numpy.typing import ArrayLike
import matplotlib.pyplot as plt
import pandas as pd
import yaml
from yaml import CLoader as Loader
import pandas as pd
from phonopy import PhonopyQHA

def qha_properties(
    ev_csv: str,
    phonopy_save_path: str,
    thermal_properties_file_pattern: str,
    pressure: Optional[float] = None,
    t_min: float = 0,
    t_step: float = 10,
    t_max: float = 1000,
) -> Dict:

    ev_df = pd.read_csv(ev_csv)
    volumes = ev_df['volumes'].to_numpy()
    energies = ev_df['energies'].to_numpy()
    thermal_properties_files = glob(thermal_properties_file_pattern)
    thermal_properties_files = sorted(thermal_properties_files)
    assert len(thermal_properties_files) > 0, \
        f'No files matching {thermal_properties_file_pattern} in {os.getcwd()}'

    entropy = []
    cv = []
    fe = []
    for filename in thermal_properties_files:
        print("Reading %s" % filename)
        with open(filename) as fp:
            tp_results = yaml.load(fp, Loader=Loader)
        thermal_properties = tp_results["thermal_properties"]
        temperatures = [v["temperature"] for v in thermal_properties]
        cv.append([v["heat_capacity"] for v in thermal_properties])
        entropy.append([v["entropy"] for v in thermal_properties])
        fe.append([v["free_energy"] for v in thermal_properties])

    qha = PhonopyQHA(
        volumes,
        energies,
        temperatures=temperatures,
        pressure=pressure,
        free_energy=np.transpose(fe),
        cv=np.transpose(cv),
        entropy=np.transpose(entropy),
        t_max=t_max,
        verbose=True,
    )

    try:
        qha.plot_helmholtz_volume()
        qha.write_helmholtz_volume()
        plt.savefig('helmholtz_volume.png')
    except:
        print('Failed to calculate helmholtz_volume')
    try:
        qha.plot_volume_temperature()
        plt.savefig('volume_temperature.png')
        qha.write_volume_temperature()
    except:
        print('Failed to calculate volume_temperature')
    try:
        qha.plot_thermal_expansion()
        plt.savefig('thermal_expansion.png')
        qha.write_thermal_expansion()
    except:
        print('Failed to calculate thermal_expansion')
    try:
        qha.plot_gibbs_temperature()
        plt.savefig('gibbs_temperature.png')
        qha.write_gibbs_temperature()
    except:
        print('Failed to calculate gibbs_temperature')
    try:
        qha.plot_bulk_modulus_temperature()
        plt.savefig('bulk_modulus_temperature.png')
        qha.write_bulk_modulus_temperature()
    except:
        print('Failed to calculate bulk_modulus_temperature')
    try:
        qha.plot_heat_capacity_P_numerical()
        plt.savefig('heat_capacity_P_numerical.png')
        qha.write_heat_capacity_P_numerical()
    except:
        print('Failed to calculate heat_capacity_P_numerical')
    try:
        qha.plot_heat_capacity_P_polyfit()
        plt.savefig('heat_capacity_P_polyfit.png')
        qha.write_heat_capacity_P_polyfit()
    except:
        print('Failed to calculate heat_capacity_P_polyfit')
    try:
        qha.plot_gruneisen_temperature()
        plt.savefig('gruneisen_temperature.png')
        qha.write_gruneisen_temperature()
    except:
        print('Failed to calculate gruneisen_temperature')

    return {}
