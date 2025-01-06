#!/usr/bin/env python
'''
Asimmodule for training MACE models

Author: mkphuthi@github.com
'''
from typing import Dict, Optional, Union
import sys
from pathlib import Path
import logging
import warnings
warnings.filterwarnings("ignore")
import json
from numpy.random import randint
from mace.cli.run_train import main as mace_run_train_main
from mace.cli.create_lammps_model import main as create_lammps_model

def train_mace(
    config: Union[Dict,str],
    randomize_seed: bool = False,
    compile_lammps: bool = False,
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
            config = json.load(fp)
    
    if randomize_seed:
        config['seed'] = randint(0, 1000000)

    config_file_path = str(Path("mace_config.yaml").resolve())
    with open(config_file_path, "w") as f:
        json.dump(config, f, indent=2)

    logging.getLogger().handlers.clear()
    sys.argv = ["program", "--config", config_file_path]
    mace_run_train_main()

    if compile_lammps:
        create_lammps_model('mace_test_compiled.model')
    return {}
