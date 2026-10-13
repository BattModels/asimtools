#!/usr/bin/env python
'''
Sets up and runs a VASP NEB (Nudged Elastic Band) calculation given an initial
band of images. Heavily uses pymatgen for IO and MP settings.
VASP must be installed.

For a VASP NEB calculation, each image of the band is placed in its own
numbered directory (``00``, ``01``, ..., ``0N``) where ``00`` and ``0N`` are
the fixed endpoints and the directories in between hold the intermediate
images. The shared INCAR, KPOINTS and POTCAR live in the root directory and
the ``IMAGES`` INCAR tag is set to the number of intermediate images.

VASP does not compute the endpoints in NEB mode, so ``00/OUTCAR`` and
``0N/OUTCAR`` must hold converged endpoint calculations. With
``relax_endpoints`` they are produced here by relaxing the endpoints with
:func:`asimtools.asimmodules.vasp.vasp.vasp`; otherwise they must already exist
(e.g. from a previous run) or an exception is raised.

Author: mkphuthi@github.com
'''
import os
import shutil
import logging
from copy import deepcopy
import numpy as np
from ase.io import read
from ase.geometry import find_mic
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp import Poscar
from asimtools.utils import get_images
from asimtools.asimmodules.vasp.utils import write_vasp_inputs
from asimtools.asimmodules.vasp.vasp import execute_vasp_run_command, vasp

NEB_ONLY_TAGS = ('IMAGES', 'ICHAIN', 'SPRING', 'LCLIMB', 'LTANGENTOLD', 'LNEBCELL')

def vasp_neb(
    images: dict,
    user_incar_settings: dict | None = None,
    user_kpoints_settings: dict | None = None,
    user_potcar_functional: str = 'PBE_64',
    potcar: dict | None = None,
    vaspinput_kwargs: dict | None = None,
    command: str = 'srun vasp_std',
    mpset: str | None = None,
    prev_calc: os.PathLike | None = None,
    image_output_file: str | None = 'OUTCAR',
    run_neb: bool = True,
    relax_endpoints: bool = False,
    relax_vasp_args: dict | None = None,
) -> dict:
    """Set up and run a VASP NEB calculation for the given initial band

    :param images: Initial band of images for the NEB calculation. The first
        and last images are treated as the fixed endpoints and the rest as the
        intermediate images. Images specification, see
        :func:`asimtools.utils.get_images`
    :type images: dict
    :param user_incar_settings: Dictionary of INCAR settings to override
        defaults or MP settings, defaults to None. The ``IMAGES`` tag is set
        automatically to the number of intermediate images unless explicitly
        provided here.
    :type user_incar_settings: dict, optional
    :param user_kpoints_settings: Dictionary of KPOINTS settings to override
        defaults or MP settings, defaults to None
    :type user_kpoints_settings: dict, optional
    :param user_potcar_functional: Potcar functional to use in case of MP
        settings, defaults to 'PBE_64'
    :type user_potcar_functional: str, optional
    :param potcar: Dictionary specifying Potcar settings, see
        :class:`pymatgen.io.vasp.inputs.Potcar`, defaults to None
    :type potcar: dict, optional
    :param vaspinput_kwargs: Dictionary of pymatgen's VaspInput arguments.
        See :class:`pymatgen.io.vasp.inputs.VaspInput`
    :type vaspinput_kwargs: dict, optional
    :param command: Command with which to run VASP, defaults to 'srun vasp_std'
    :type command: str, optional
    :param mpset: Materials Project VASP set to use see
        :mod:`pymatgen.io.vasp.sets`, defaults to None
    :type mpset: str, optional
    :param prev_calc: Path to previous VASP calculation to use as starting
        point for MP settings, defaults to None
    :type prev_calc: os.PathLike, optional
    :param image_output_file: File to read output image from in each image
        directory and write in standard asimtools format, set to None to skip,
        defaults to 'OUTCAR'
    :type image_output_file: str, optional
    :param run_neb: Whether to run the VASP NEB after writing input files,
        defaults to True
    :type run_neb: bool, optional
    :param relax_endpoints: Relax the first and last images with
        :func:`asimtools.asimmodules.vasp.vasp.vasp` in ``relax_initial`` and
        ``relax_final`` before setting up the NEB, defaults to False. If False,
        converged endpoint OUTCARs must already exist in the endpoint
        directories
    :type relax_endpoints: bool, optional
    :param relax_vasp_args: Arguments passed to
        :func:`asimtools.asimmodules.vasp.vasp.vasp` for the endpoint
        relaxations. By default the relaxations use the same arguments as the
        NEB with the NEB-only INCAR tags removed; ``user_incar_settings`` given
        here is merged into those INCAR settings and any other key replaces
        the NEB value. The endpoint convergence check uses ``-EDIFFG`` of the
        resulting relaxation INCAR settings, defaults to None
    :type relax_vasp_args: dict, optional
    """

    band = get_images(**images)
    assert len(band) >= 3, \
        'NEB requires at least 3 images (2 endpoints + 1 intermediate)'

    n_intermediate = len(band) - 2
    image_dirs = [f'{i:02d}' for i in range(len(band))]
    user_incar_settings = deepcopy(user_incar_settings) or {}

    relax_vasp_args = deepcopy(relax_vasp_args) or {}
    for key in ('image', 'run_vasp'):
        if key in relax_vasp_args:
            raise ValueError(f'{key} is set by vasp_neb, not relax_vasp_args')
    relax_incar = {
        k: v for k, v in user_incar_settings.items() if k not in NEB_ONLY_TAGS
    }
    relax_incar.update(relax_vasp_args.pop('user_incar_settings', None) or {})
    relax_args = {
        'user_incar_settings': relax_incar,
        'user_kpoints_settings': user_kpoints_settings,
        'user_potcar_functional': user_potcar_functional,
        'potcar': potcar,
        'vaspinput_kwargs': vaspinput_kwargs,
        'command': command,
        'mpset': mpset,
        'prev_calc': prev_calc,
        'image_output_file': image_output_file,
    }
    relax_args.update(relax_vasp_args)

    ediffg = relax_incar.get('EDIFFG')
    if ediffg is None or ediffg >= 0:
        raise ValueError(
            'Endpoint convergence check needs a force criterion: set a '
            f'negative EDIFFG for the relaxation, got {ediffg}'
        )
    fmax_tol = -ediffg

    endpoints = ((0, 'relax_initial'), (-1, 'relax_final'))
    if relax_endpoints:
        for idx, relax_dir in endpoints:
            os.makedirs(relax_dir, exist_ok=True)
            cwd = os.getcwd()
            os.chdir(relax_dir)
            try:
                vasp({'atoms': band[idx]}, run_vasp=True, **relax_args)
            finally:
                os.chdir(cwd)
            os.makedirs(image_dirs[idx], exist_ok=True)
            shutil.copy(
                os.path.join(relax_dir, 'OUTCAR'),
                os.path.join(image_dirs[idx], 'OUTCAR'),
            )

    # Endpoints must have converged OUTCARs whether relaxed now or previously
    for idx, _ in endpoints:
        outcar = os.path.join(image_dirs[idx], 'OUTCAR')
        if not os.path.exists(outcar):
            raise FileNotFoundError(
                f'No endpoint OUTCAR found at {outcar}, set relax_endpoints '
                'to True to compute it'
            )
        relaxed = read(outcar)
        relaxed.set_constraint(band[idx].constraints)
        fmax = np.linalg.norm(relaxed.get_forces(), axis=1).max()
        if fmax > fmax_tol:
            raise RuntimeError(
                f'Endpoint {outcar} is not converged: fmax {fmax:.4f} > '
                f'{fmax_tol} eV/Ang'
            )
        # Undo VASP's wrapping so the endpoint stays continuous with the band
        disp, _ = find_mic(
            relaxed.positions - band[idx].positions, band[idx].cell, pbc=True
        )
        relaxed.positions = band[idx].positions + disp
        band[idx] = relaxed

    # Set the IMAGES INCAR tag to the number of intermediate images unless the
    # user explicitly specified it
    user_incar_settings.setdefault('IMAGES', n_intermediate)

    # Write the shared INCAR, KPOINTS and POTCAR to the root directory using
    # the first endpoint as the reference image
    write_vasp_inputs(
        {'atoms': band[0]},
        mpset=mpset,
        user_incar_settings=user_incar_settings,
        user_kpoints_settings=user_kpoints_settings,
        user_potcar_functional=user_potcar_functional,
        potcar=potcar,
        vaspinput_kwargs=vaspinput_kwargs,
        prev_calc=prev_calc,
    )

    # VASP NEB does not use the root POSCAR, instead each image gets its own
    # numbered directory with a POSCAR
    if os.path.exists('POSCAR'):
        os.remove('POSCAR')

    for image_dir, atoms in zip(image_dirs, band):
        os.makedirs(image_dir, exist_ok=True)
        struct = AseAtomsAdaptor.get_structure(atoms)
        Poscar(struct).write_file(os.path.join(image_dir, 'POSCAR'))

    if run_neb:
        execute_vasp_run_command(command)

        if image_output_file:
            for image_dir in image_dirs:
                output_path = os.path.join(image_dir, image_output_file)
                if os.path.exists(output_path):
                    atoms_output = read(output_path)
                    atoms_output.write(
                        os.path.join(image_dir, 'image_output.xyz'),
                        format='extxyz',
                    )
                else:
                    logging.warning(
                        'No %s found in %s, skipping output for this image',
                        image_output_file, image_dir,
                    )

    return {}
