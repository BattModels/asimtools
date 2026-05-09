#!/usr/bin/env python
'''
Runs VASP based on input files and optionally MP settings. 
Heavily uses pymatgen for IO and MP settings.
VASP must be installed 

Author: mkphuthi@github.com
'''
from typing import Dict, Optional
import os
import logging
from ase.io import read
from pymatgen.io.vasp import Poscar, Incar, Potcar, Kpoints, VaspInput
import pymatgen.io.vasp.sets
from custodian.custodian import Custodian
import custodian.vasp.handlers
from custodian.vasp.jobs import VaspJob
from asimtools.utils import (
    get_atoms,
)

def execute_vasp_run_command(
    command: str,
    max_errors: int = 10,
    vasp_job_kwargs: Optional[Dict] = None,
    handlers: list[str] = ["VaspErrorHandler", "UnconvergedErrorHandler"],
    handler_kwargs: list[dict] = None,
) -> None:
    command = command.split(' ')
    if vasp_job_kwargs is None:
        vasp_job_kwargs = {}

    vaspjob = VaspJob(
        vasp_cmd=command,
        output_file="vasp_stdout.txt",
        stderr_file="vasp_stderr.txt",
        final=vasp_job_kwargs.pop('final', True),
        backup=vasp_job_kwargs.pop('backup', False),
        **vasp_job_kwargs
    )

    if handler_kwargs is not None:
        assert len(handlers) == len(handler_kwargs), \
            "Num. handlers ({len(handlers)}) must equal num. handler_kwargs " \
            "({len(handler_kwargs)})"
    else:
        handler_kwargs = [{} for handler in handlers]

    handler_objs = []
    for i, handler in enumerate(handlers):
        try:
            handler_cls = getattr(custodian.vasp.handlers, handler)
        except AttributeError as e:
            print(f"{handler} isn't valid handler from ustodian.vasp.handlers")
        handler_objs.append(handler_cls(
            output_filename="vasp_stdout.txt",
            **handler_kwargs[i]
        ))

    jobs = [vaspjob]
    c = Custodian(handler_objs, jobs, max_errors=max_errors)
    custodian_output = c.run()

    if custodian_output[0]['nonzero_return_code']:
        err_txt = f'VASP failed with custodian'
        err_txt += '\nSee vasp_stderr.txt, OUTCAR etc. for details.'
        logging.error(err_txt)
        with open('vasp_stderr.txt', 'a+', encoding='utf-8') as f:
            f.write(str(custodian_output))

    else:
        logging.info(
            f'VASP run completed.'
        )

def custodian_vasp(
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
    max_errors: int = 10,
    vasp_job_kwargs: Optional[Dict] = None,
    handlers: list[str] = ["VaspErrorHandler", "UnconvergedErrorHandler"]
) -> Dict:
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
    :param write_image_output: Whether to write output image in standard 
        asimtools format to file, defaults to True
    :type write_image_output: bool, optional
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

    if run_vasp:
        execute_vasp_run_command(
            command,
            max_errors=max_errors,
            vasp_job_kwargs=vasp_job_kwargs,
            handlers=handlers,
        )

        if write_image_output:
            atoms_output = read('OUTCAR')
            atoms_output.write(
                'image_output.xyz',
                format='extxyz',
            )

    return {}
