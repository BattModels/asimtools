'''.Xauthority'''

from typing import Dict, Tuple, List
import numpy as np
import matplotlib.pyplot as plt
from ase.eos import EquationOfState
from asimtools.utils import (
    get_atoms
)

@step(job_type='process')
def preprocess(
    image: Dict,
    scale_range: Tuple[float, float] = (0.95,1.05),
    nimages: int = 5
) -> Tuple[bool,List,List]:
    ''' Get the images for the eos calculation '''
    # Fail early principle. There should usually be assert statements here
    # to make sure inputs are correct
    assert scale_range[0] < scale_range[1], 'x_scale[0] should be \
        smaller than x_scale[1]'

    atoms = get_atoms(**image)

    # Preprocessing the data, structures and data structures
    scales = np.linspace(scale_range[0], scale_range[1], nimages)
    images = []
    for scale in scales:
        new_cell = atoms.get_cell() * scale
        new_atoms = atoms.copy()
        new_atoms.set_cell(new_cell)
        images.append(new_atoms)

    return True, singlepoint_sim_inputs

@step(job_type='simulation')
def calculate_energies(images, sim_input, scales):
    ''' Performs the energy calculations '''
    distribute_script(images, sim_inputs)

@step(job_type='process')
def postprocess(eos_output_csv):
    ''' plot things '''
    data = np.genfromtxt(eos_output_csv, delimiter=',')
    volumes = data[:,1]
    energies = data[:,0]
    eos_fit = EquationOfState(volumes, energies)
    try:
        v0, e0, B = eos_fit.fit()
    except ValueError:
        print('ERROR: Could not fit EOS, check structures')
        v0, e0, B = None, None, None
        job.update_status('failed')

    fig, ax = plt.subplots()
    if v0 is not None:
        B_GPa = B / kJ * 1.0e24 # eV/Ang^3 to GPa
        eos_fit.plot(ax=ax)

        # Get equilibrium lattice scaling by interpolation
        xs2 = np.linspace(scale_range[0], scale_range[1], 100)
        vs2 = []
        for i in xs2:
            samp_atoms = atoms.copy()
            samp_atoms.cell = samp_atoms.cell*i
            vs2.append(samp_atoms.get_volume())
        x0 = np.interp(v0, vs2, xs2)

        job.update_output({
            'equilibrium_volume': float(v0),
            'equilibrium_energy': float(e0),
            'equilibrium_scale': float(x0),
            'bulk_modulus': float(B_GPa),
        })
    else:
        ax.plot(volumes, energies)
    ax.set_xlabel(r'Volume ($\AA$)')
    ax.set_ylabel(r'Energy (eV)')
    plt.savefig(join_names([job.get_prefix(), 'eos.png']))
    plt.close(fig)

def eos(
    calc_input: Dict,
    image: dict,
    prefix: str = '',
    nimages: int = 5,
    scale_range: Tuple[float] = (0.95, 1.05),
    workdir: str = '.',
    **kwargs
) -> Dict:

    status, images, scales = preprocess(image, scale_range, nimages)
    sim_results = calculate_energies(images, scales)
    analysis_results = postprocess(sim_results)

