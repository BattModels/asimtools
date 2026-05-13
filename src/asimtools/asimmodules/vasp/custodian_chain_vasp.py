#!/usr/bin/env python
'''
Runs VASP based on input files and optionally MP settings. 
Heavily uses pymatgen for IO and MP settings.
VASP must be installed 

Author: mkphuthi@github.com
'''
import os
import logging
from ase.io import read
from custodian.custodian import Custodian
import custodian.vasp.handlers
from custodian.vasp.jobs import VaspJob
from asimtools.asimmodules.vasp.utils import write_vasp_inputs

def execute_vasp_run_command(
    command: str,
    max_errors: int = 10,
    vasp_job_kwargs_set: list[dict] | None = None,
    handlers: list[str] | None = None,
    handler_kwargs: list[dict] | None = None,
) -> None:
    command = command.split(' ')
    if handlers is None:
        handlers = ["VaspErrorHandler", "UnconvergedErrorHandler"]

    jobs = []
    for i, vasp_job_kwargs in enumerate(vasp_job_kwargs_set):
        final = False
        if i == len(vasp_job_kwargs_set) - 1:
            final = True

        vaspjob = VaspJob(
            vasp_cmd=command,
            output_file="vasp_stdout.txt",
            stderr_file="vasp_stderr.txt",
            final=vasp_job_kwargs.get('final', final),
            backup=vasp_job_kwargs.get('backup', False),
            **{k: v for k, v in vasp_job_kwargs.items() if k not in ('final', 'backup')}
        )
        jobs.append(vaspjob)

    if handler_kwargs is not None:
        if len(handlers) != len(handler_kwargs):
            raise ValueError(
                f"Num. handlers ({len(handlers)}) must equal num. "
                f"handler_kwargs ({len(handler_kwargs)})"
            )
    else:
        handler_kwargs = [{} for handler in handlers]

    handler_objs = []
    for i, handler in enumerate(handlers):
        try:
            handler_cls = getattr(custodian.vasp.handlers, handler)
        except AttributeError:
            raise ValueError(f"{handler} is not a valid handler from custodian.vasp.handlers")
        handler_objs.append(handler_cls(
            output_filename="vasp_stdout.txt",
            **handler_kwargs[i]
        ))

    c = Custodian(handler_objs, jobs, max_errors=max_errors)
    custodian_output = c.run()

    if any(out['nonzero_return_code'] for out in custodian_output):
        err_txt = f'VASP failed with custodian'
        err_txt += '\nSee vasp_stderr.txt, OUTCAR etc. for details.'
        logging.error(err_txt)
        with open('vasp_stderr.txt', 'a+', encoding='utf-8') as f:
            f.write(str(custodian_output))

    else:
        logging.info(
            f'VASP run completed.'
        )

def custodian_chain_vasp(
    image: dict | None,
    vasp_job_kwargs_set: list[dict],
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
    max_errors: int = 10,
    handlers: list[str] | None = None,
    handler_kwargs: list[dict] | None = None,
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
    :param handler_kwargs: List of kwargs dicts for each handler, must match
        length of handlers, defaults to None
    :type handler_kwargs: list[dict], optional
    :param run_vasp: Whether to run VASP after writing input files,
        defaults to True
    :type run_vasp: bool, optional
    """

    if len(vasp_job_kwargs_set) == 0:
        raise ValueError('Must supply vasp_job_kwargs_set')

    write_vasp_inputs(
        image,
        mpset=mpset,
        user_incar_settings=user_incar_settings,
        user_kpoints_settings=user_kpoints_settings,
        user_potcar_functional=user_potcar_functional,
        potcar=potcar,
        vaspinput_kwargs=vaspinput_kwargs,
        prev_calc=prev_calc,
    )

    if run_vasp:
        execute_vasp_run_command(
            command,
            max_errors=max_errors,
            vasp_job_kwargs_set=vasp_job_kwargs_set,
            handlers=handlers,
            handler_kwargs=handler_kwargs,
        )

        if image_output_file:
            atoms_output = read(image_output_file)
            atoms_output.write(
                'image_output.xyz',
                format='extxyz',
            )

    return {}
