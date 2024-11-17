#!/usr/bin/env python
'''
Runs a user defined lammps script or template. LAMMPS must be installed 

Author: mkphuthi@github.com
'''
from typing import Dict, Optional
import sys
from pathlib import Path
import logging
import warnings
warnings.filterwarnings("ignore")
from numpy.random import randint
from mace.cli.run_train import main as mace_run_train_main

def train_mace(
    config: Union[Dict,str],
    randomize_seed: bool = False,
) -> Dict:
    """Runs MACE training

    :param config: MACE config dictionary or path to config file
    :type config: Union[Dict,str]
    :param randomize_seed: Whether to randomize the seed, defaults to False
    :type randomize_seed: bool
    :return: Dictionary of results
    :rtype: Dict
    """

    if isinstance(config, str):
        with open(config, 'r') as fp:
            config = yaml.safe_load(fp)
    
    if randomize_seed:
        config['seed'] = randint(0, 1000000)

    config_file_path = str(Path("mace_config.yaml").resolve())
    with open(config_file_path, "w") as f:
        f.write(config)

    logging.getLogger().handlers.clear()
    sys.argv = ["program", "--config", config_file_path]
    mace_run_train_main()
    return {}


train_mace(
    config={
        "model": "MACE",
        "num_channels": 32,
        "max_L": 0,
        "r_max": 4.0,
        "name": "mace02_com1_gen1",
        "model_dir": "MACE_models",
        "log_dir": "MACE_models",
        "checkpoints_dir": "MACE_models",
        "results_dir": "MACE_models",
        "train_file": "data/solvent_xtb_train_23_gen1.xyz",
        "valid_file": "data/solvent_xtb_train_50.xyz",
        "test_file": "data/solvent_xtb_test.xyz",
        "energy_key": "energy_xtb",
        "forces_key": "forces_xtb",
        "device": "cuda",
        "batch_size": 10,
        "max_num_epochs": 500,
        "swa": True,
        "seed": 1234,
    }
)
