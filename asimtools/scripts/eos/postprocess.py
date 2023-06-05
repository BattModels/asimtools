'''.Xauthority'''

from typing import Dict, Tuple
import numpy as np
import matplotlib.pyplot as plt
from ase.eos import EquationOfState
from ase. units import kJ
from asimtools.job import leaf, load_output_images, load_input_images
from asimtools.utils import (
    write_csv_from_dict
)

@leaf
def postprocess(
    config_input: Dict,
    scale_range: Tuple[float, float],
    **kwargs
) -> Tuple[None,Dict]:
    ''' plot things '''
    # Should we standardize to using pandas? Issue is that switching
    # between np and pandas is probably not good

    # Should change this to load_jobs_from_directory once we
    # know everything else works
    images = load_output_images(pattern='../step-0/id-*')
    volumes = np.array([at.get_volume() for at in images])
    energies = np.array([at.get_potential_energy() for at in images])
    write_csv_from_dict(
        'eos_output.csv',
        {'volumes': volumes, 'energies': energies}
    )
    eos_fit = EquationOfState(volumes, energies)
    print('here')
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
        atoms = load_input_images('../step-0')[0]
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
    return None, results
