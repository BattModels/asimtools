#!/usr/bin/env python
'''
Runs VASP based on input files and optionally MP settings. 
Heavily uses pymatgen for IO and MP settings.
VASP must be installed 

Author: mkphuthi@github.com
'''
from typing import Dict, Optional, Sequence
import os
import sys
from pathlib import Path
from numpy.random import randint
import subprocess
import logging
from ase.io import read
from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp import Poscar, Incar, Potcar, Kpoints, VaspInput
import pymatgen.io.vasp.sets
from asimtools.utils import (
    get_atoms,
)

def vasp(
    image: Optional[Dict],
    user_incar_settings: Optional[Dict] = None,
    user_kpoints_settings: Optional[Dict] = None,
    user_potcar_functional: str = 'PBE_64',
    potcar: Optional[Dict] = None,
    vaspinput_kwargs: Optional[Dict] = None,
    command: str = 'srun vasp_std',
    mpset: Optional[str] = None,
    prev_calc: Optional[os.PathLike] = None,
    write_image_output: bool = True,
    run_vasp: bool = True,
) -> Dict:
    """Run VASP with given input files and specified image

    :param image: Initial image for VASP calculation. Image specification,
        see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param vaspinput_args: Dictionary of pymatgen's VaspInput arguments. 
        See :class:`pymatgen.io.vasp.inputs.VaspInput`
    :type vaspinput_args: Dict
    :param command: Command with which to run VASP, defaults to 'vasp_std'
    :type command: str, optional
    :param mpset: Materials Project VASP set to use see 
        :mod:`pymatgen.io.vasp.sets`, defaults to None
    :type mpset: str, optional
    :param write_image_output: Whether to write output image in standard 
        asimtools format to file, defaults to False
    :type write_image_output: bool, optional
    """

    if vaspinput_kwargs is None:
        vaspinput_kwargs = {}

    struct = get_atoms(**image, return_type='pymatgen')
    if mpset is not None:
        # if not ((incar is None) or (potcar is None) or (kpoints is None)):
        #     raise ValueError(
        #         'Provide either mpset or all of incar and kpoints'
        #     )
        try:
            set_ = getattr(pymatgen.io.vasp.sets, mpset)
        except:
            raise ImportError(
                f'Unknown mpset: {mpset}. See available sets in pymatgen.')

        if prev_calc is not None:
            vasp_input = set_.from_prev_calc(
                prev_calc,
                user_incar_settings=user_incar_settings,
                user_kpoints_settings=user_kpoints_settings,
                user_potcar_functional=user_potcar_functional,
                **vaspinput_kwargs
            )
        else:
            vasp_input = set_(
                struct,
                user_incar_settings=user_incar_settings,
                user_kpoints_settings=user_kpoints_settings,
                user_potcar_functional=user_potcar_functional,
                **vaspinput_kwargs
            )
        
    else:

        incar = Incar(user_incar_settings)
        incar.check_params()
        if potcar is not None:
            potcar = Potcar(potcar)
        if user_kpoints_settings is not None:
            kpoints = Kpoints(user_kpoints_settings)
        else:
            kpoints=None

        if vaspinput_args is None:
            vaspinput_args = {}

        vasp_input = VaspInput(
            incar=incar,
            kpoints=kpoints,
            poscar=Poscar(struct),
            potcar=potcar,
            **vaspinput_kwargs
        )


    vasp_input.write_input("./")
    # if incar_kwargs is not None:
    #     with open('INCAR', 'a+') as fp:
    #         for k, v in incar_kwargs.items():
    #             fp.write(f'\n{k} = {v}')

    if run_vasp:
        command = command.split(' ')
        completed_process = subprocess.run(
                command, check=False, capture_output=True, text=True,
            )

        with open('vasp_stdout.txt', 'a+', encoding='utf-8') as f:
            f.write(completed_process.stdout)

        if completed_process.returncode != 0:
            err_txt = f'Failed to run VASP\n'
            err_txt += 'See vasp_stderr.txt for details.'
            logging.error(err_txt)
            with open('vasp_stderr.txt', 'a+', encoding='utf-8') as f:
                f.write(completed_process.stderr)
            completed_process.check_returncode()
            return {}

        if write_image_output:
            atoms_output = read('OUTCAR')
            atoms_output.write(
                'image_output.xyz',
                format='extxyz',
            )

    return {}
