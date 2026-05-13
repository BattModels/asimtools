#!/usr/bin/env python
'''
Runs VASP based on input files and optionally MP settings. 
Heavily uses pymatgen for IO and MP settings.
VASP must be installed 

Author: mkphuthi@github.com
'''
import os
import subprocess
import logging
import shutil
from ase.io import read
from pymatgen.io.vasp import Poscar, Incar, Potcar, Kpoints, VaspInput
import pymatgen.io.vasp.sets
from asimtools.utils import (
    get_atoms,
    get_str_btn,
)

def execute_vasp_run_command(
    command: str,
) -> None:
    command = command.split(' ')
    completed_process = subprocess.run(
            command, check=False, capture_output=True, text=True,
        )

    with open('vasp_stdout.txt', 'a+', encoding='utf-8') as f:
        f.write(completed_process.stdout)

    if completed_process.returncode != 0:
        err_txt = f'VASP failed with error code: {completed_process.returncode}'
        err_txt += '\nSee vasp_stderr.txt for details.'
        logging.error(err_txt)
        with open('vasp_stderr.txt', 'a+', encoding='utf-8') as f:
            f.write(completed_process.stderr)
        completed_process.check_returncode()
        if "ZBRENT: fatal error in bracketing" in completed_process.stderr:
            filenum = 0
            while os.path.exists(f'OUTCAR.{filenum}'):
                filenum += 1
            shutil.move('OUTCAR', f'OUTCAR.{filenum}')
            shutil.move('POSCAR', f'POSCAR.{filenum}')
            shutil.move('CONTCAR', f'POSCAR')
            execute_vasp_run_command(command)
    else:
        logging.info(
            f'VASP run completed with code {completed_process.returncode}.'
        )

def vasp(
    image: dict | None,
    user_incar_settings: dict | None = None,
    user_kpoints_settings: dict | None = None,
    user_potcar_functional: str = 'PBE_64',
    potcar: dict | None = None,
    vaspinput_kwargs: dict | None = None,
    command: str = 'srun vasp_std',
    mpset: str | None = None,
    prev_calc: os.PathLike | None = None,
    image_output_file: str | None = 'OUTCAR',
    run_vasp: bool = True,
) -> dict:
    """Run VASP with given input files and specified image

    :param image: Initial image for VASP calculation. Image specification,
        see :func:`asimtools.utils.get_atoms`
    :type image: Dict
    :param user_incar_settings: Dictionary of INCAR settings to override
        defaults or MP settings, defaults to None
    :type user_incar_settings: Dict, optional
    :param user_kpoints_settings: Dictionary of KPOINTS settings to override
        defaults or MP settings, defaults to None
    :type user_kpoints_settings: Dict, optional
    :param user_potcar_functional: Potcar functional to use in case of MP
        settings, defaults to 'PBE_64'
    :type user_potcar_functional: str, optional
    :param potcar: Dictionary specifying Potcar settings, see
        :class:`pymatgen.io.vasp.inputs.Potcar`, defaults to None
    :type potcar: Dict, optional
    :param prev_calc: Path to previous VASP calculation to use as starting
        point for MP settings, defaults to None
    :type prev_calc: os.PathLike, optional
    :param vaspinput_kwargs: Dictionary of pymatgen's VaspInput arguments. 
        See :class:`pymatgen.io.vasp.inputs.VaspInput`
    :type vaspinput_kwargs: Dict
    :param command: Command with which to run VASP, defaults to 'vasp_std'
    :type command: str, optional
    :param mpset: Materials Project VASP set to use see 
        :mod:`pymatgen.io.vasp.sets`, defaults to None
    :type mpset: str, optional
    :param image_output_file: File to read output image from and write in
        standard asimtools format, set to None to skip, defaults to 'OUTCAR'
    :type image_output_file: str, optional
    :param run_vasp: Whether to run VASP after writing input files,
        defaults to True
    :type run_vasp: bool, optional
    """

    if vaspinput_kwargs is None:
        vaspinput_kwargs = {}

    struct = get_atoms(**image, return_type='pymatgen')
    if mpset is not None:
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

        vasp_input = VaspInput(
            incar=incar,
            kpoints=kpoints,
            poscar=Poscar(struct),
            potcar=potcar,
            **vaspinput_kwargs
        )

    vasp_input.write_input("./")

    if run_vasp:
        optimization_failed = False
        execute_vasp_run_command(command)
        incar = vasp_input.incar
        ibrion = incar.get('IBRION', -1)
        nsw = incar.get('NSW', 0)
        # Don't write result if running a relaxation that didn't finish
        if ibrion in (1,2,3):
            if nsw > 0:
                with open('OUTCAR', 'r') as f:
                    lines = f.readlines()
                lines = lines[::-1]
                for line in lines:
                    if 'Ionic step' in line:
                        last_step = int(get_str_btn(line, 'Ionic step', '--'))
                        if last_step == nsw:
                            optimization_failed = True
                            raise RuntimeError(
                                'VASP relaxation did not complete. '
                                'Check OUTCAR for details.'
                            )

        if image_output_file and not optimization_failed:
            atoms_output = read(image_output_file)
            atoms_output.write(
                'image_output.xyz',
                format='extxyz',
            )

    return {}
