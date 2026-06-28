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

Author: mkphuthi@github.com
'''
import os
import logging
from copy import deepcopy
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp import Poscar
from asimtools.utils import get_images
from asimtools.asimmodules.vasp.utils import write_vasp_inputs
from asimtools.asimmodules.vasp.vasp import execute_vasp_run_command

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
    run_vasp: bool = True,
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
    :param run_vasp: Whether to run VASP after writing input files,
        defaults to True
    :type run_vasp: bool, optional
    """

    band = get_images(**images)
    assert len(band) >= 3, \
        'NEB requires at least 3 images (2 endpoints + 1 intermediate)'

    n_intermediate = len(band) - 2

    # Set the IMAGES INCAR tag to the number of intermediate images unless the
    # user explicitly specified it
    user_incar_settings = deepcopy(user_incar_settings) or {}
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

    image_dirs = []
    for i, atoms in enumerate(band):
        image_dir = f'{i:02d}'
        os.makedirs(image_dir, exist_ok=True)
        struct = AseAtomsAdaptor.get_structure(atoms)
        Poscar(struct).write_file(os.path.join(image_dir, 'POSCAR'))
        image_dirs.append(image_dir)

    if run_vasp:
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
