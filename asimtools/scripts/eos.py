'''.Xauthority'''

from typing import Dict, Tuple, List
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
from ase.eos import EquationOfState
from ase. units import kJ
from asimtools.job import UnitJob
from asimtools.utils import (
    get_atoms,
)

def singlepoint_calculations(
    config_input: Dict,
    image: Dict,
    image_array_input: Dict,
    nimages: int = 5,
    scale_range: Tuple[float] = (0.95, 1.05),
):
    # Fail early principle. There should usually be assert statements 
    # here to make sure inputs are correct
    assert scale_range[0] < scale_range[1], 'x_scale[0] should be smaller than\
         x_scale[1]'

    atoms = get_atoms(**image)

    # Preprocessing the data, structures and data structures
    scales = np.linspace(scale_range[0], scale_range[1], nimages)

    # single point calculation script called here. Note that this will
    # automatically submit nimages jobs
    scaled_atoms = []
    scaled_ids = []
    for scale in scales:
        new_cell = atoms.get_cell() * scale
        new_atoms = atoms.copy()
        new_atoms.set_cell(new_cell)
        scaled_atoms.append(new_atoms)
        scaled_ids.append(f'x{scale:.2f}')

    image_array_input = deepcopy(image_array_input)
    image_array_job = UnitJob(config_input, image_array_input)
    image_array_job.submit()

def eos_postprocess(
    config_input: Dict,
    eos_output_csv: str,
    scale_range: Tuple[float, float],
    eos_output_image_file: str, 
):
    ''' plot things '''
    # Should we standardize to using pandas? Issue is that switching
    # between np and pandas is probably not good
    data = np.genfromtxt(eos_output_csv, delimiter=',')
    volumes = data[:,1]
    energies = data[:,0]
    eos_fit = EquationOfState(volumes, energies)
    try:
        v0, e0, B = eos_fit.fit()
    except ValueError:
        print('ERROR: Could not fit EOS, check structures')
        v0, e0, B = None, None, None

    fig, ax = plt.subplots()
    if v0 is not None:
        B_GPa = B / kJ * 1.0e24 # eV/Ang^3 to GPa
        eos_fit.plot(ax=ax)

        # Get equilibrium lattice scaling by interpolation
        xs2 = np.linspace(scale_range[0], scale_range[1], 100)
        vs2 = []
        atoms = get_atoms(image_file=eos_output_image_file)
        for i in xs2:
            samp_atoms = atoms.copy()
            samp_atoms.cell = samp_atoms.cell*i
            vs2.append(samp_atoms.get_volume())
        x0 = np.interp(v0, vs2, xs2)

        results = {
            'equilibrium_volume': float(v0),
            'equilibrium_energy': float(e0),
            'equilibrium_scale': float(x0),
            'bulk_modulus': float(B_GPa),
        }
    else:
        ax.plot(volumes, energies)
        results = {}
    ax.set_xlabel(r'Volume ($\AA$)')
    ax.set_ylabel(r'Energy (eV)')
    plt.savefig('eos.png')
    plt.close(fig)
    return results

def eos(
    calc_input: Dict,
    image: dict,
    nimages: int = 5,
    scale_range: Tuple[float] = (0.95, 1.05),
    **kwargs
) -> Dict:

    chain = ChainedJobs(

    )

    chain.submit()
