#!/usr/bin/env python
'''
Shared utilities for VASP asimmodules.

Author: mkphuthi@github.com
'''
import os
from pymatgen.io.vasp import Poscar, Incar, Potcar, Kpoints, VaspInput
import pymatgen.io.vasp.sets
from asimtools.utils import get_atoms


def write_vasp_inputs(
    image: dict | None,
    mpset: str | None = None,
    user_incar_settings: dict | None = None,
    user_kpoints_settings: dict | None = None,
    user_potcar_functional: str = 'PBE_64',
    potcar: dict | None = None,
    vaspinput_kwargs: dict | None = None,
    prev_calc: os.PathLike | None = None,
) -> VaspInput:
    """Build and write VASP input files to the current directory.

    :param image: Initial image for VASP calculation. Image specification,
        see :func:`asimtools.utils.get_atoms`
    :type image: dict | None
    :param mpset: Materials Project VASP set, see
        :mod:`pymatgen.io.vasp.sets`, defaults to None
    :type mpset: str, optional
    :param user_incar_settings: INCAR settings to override defaults or MP
        settings, defaults to None
    :type user_incar_settings: dict, optional
    :param user_kpoints_settings: KPOINTS settings to override defaults or
        MP settings, defaults to None
    :type user_kpoints_settings: dict, optional
    :param user_potcar_functional: Potcar functional for MP settings,
        defaults to 'PBE_64'
    :type user_potcar_functional: str, optional
    :param potcar: Potcar settings, see
        :class:`pymatgen.io.vasp.inputs.Potcar`, defaults to None
    :type potcar: dict, optional
    :param vaspinput_kwargs: Extra kwargs for
        :class:`pymatgen.io.vasp.inputs.VaspInput`, defaults to None
    :type vaspinput_kwargs: dict, optional
    :param prev_calc: Path to previous VASP calculation for MP set
        initialisation, defaults to None
    :type prev_calc: os.PathLike, optional
    :return: The VaspInput object after writing input files
    :rtype: VaspInput
    """
    if vaspinput_kwargs is None:
        vaspinput_kwargs = {}

    struct = get_atoms(**image, return_type='pymatgen')

    if mpset is not None:
        try:
            set_ = getattr(pymatgen.io.vasp.sets, mpset)
        except AttributeError:
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
            kpoints = None

        vasp_input = VaspInput(
            incar=incar,
            kpoints=kpoints,
            poscar=Poscar(struct),
            potcar=potcar,
            **vaspinput_kwargs
        )

    vasp_input.write_input("./")
    return vasp_input
