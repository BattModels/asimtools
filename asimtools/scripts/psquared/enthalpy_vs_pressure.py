'''.Xauthority'''

from typing import Dict, Tuple
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from pymatgen.analysis.eos import EOS
from asimtools.job import branch
from asimtools.utils import get_atoms


def vinet_pv_eos(E0: float, B0: float, B1: float, V0: float):
    ''' Pressure as a function of volume based on Vinet EOS '''
    def _func(V):
        return 3 * B0 * (V/V0)**(-2/3) * (1-(V/V0)**(1/3)) * \
            np.exp(-1.5 * (B1-1) * ((V/V0)**(1/3)-1))

    return _func

@branch
def enthalpy_vs_pressure(
    config_input: Dict,
    image: Dict,
    ev_csv: str,
    plogspace: Tuple[float,float,float] = (-5,-2,10),
    sign: int = 1,
    plinspace: Tuple[float,float,float] = None,
    eos_name: str = 'vinet',
    **kwargs
) -> Tuple[None,Dict]:
    ''' Find and plot the enthalpy versus pressure given E-V data for a solid  '''

    # Get info about structure
    atoms = get_atoms(**image)
    formula = atoms.get_chemical_formula()
    try:
        spacegroup = atoms.info['spacegroup']
    except KeyError:
        spacegroup = None

    # Get E(V) data
    ev_data = np.genfromtxt(ev_csv, delimiter=',', skip_header=1)
    volumes = ev_data[:,0]
    energies = ev_data[:,1]

    # Fit EOS parameters
    eos = EOS(eos_name=eos_name).fit(volumes=volumes, energies=energies)
    ev_func = eos.func
    _, _, _, v0 = eos.eos_params
    if eos_name == 'vinet':
        pv_func = vinet_pv_eos(*eos.eos_params)

    if plinspace is not None:
        ps = sign * np.linspace(*plinspace)
    else:
        ps = sign * np.logspace(*plogspace)

    # Maybe I'll make this smarter later but this should work for now
    vmin = 0.1 * v0
    vmax = 1.2 * v0

    # Fit the chosen EOS in the spcified p-range and get h-p function
    v_samples = np.linspace(vmin, vmax, 100)
    p_samples = pv_func(v_samples)
    pv_data = np.transpose(np.vstack([v_samples, p_samples]))
    np.savetxt(
        'pressure_vs_volume.csv',
        pv_data,
        header=f'Formula:{formula}, Spacegroup:{spacegroup}\nP(eV/Ang), H(eV)',
        delimiter=','
    )

    # Plot p(V)
    _, ax = plt.subplots()
    ax.plot(v_samples, p_samples, 'o-')
    if plinspace is None:
        plt.xscale('log')
    ax.set_ylabel('Pressure (eV/$\AA^3$)')
    ax.set_xlabel('Volume ($\AA^3$)')
    plt.savefig('pressure_vs_volume.png')

    # H = U(V) + P(V)V
    h_samples = ev_func(v_samples) + pv_func(v_samples) * v_samples
    hp_func = interpolate.interp1d(p_samples, h_samples, kind='cubic')
    hp_data = np.transpose(np.vstack([ps, hp_func(ps)]))

    np.savetxt(
        'enthalpy_vs_pressure.csv',
        hp_data,
        header=f'Formula:{formula}, Spacegroup:{spacegroup}\nP(eV/Ang), H(eV)',
        delimiter=','
    )

    # Plot H(V)
    _, ax = plt.subplots()
    ax.plot(sign*hp_data[:,0], hp_data[:,1], 'o--')
    if plinspace is None:
        plt.xscale('log')
    ax.set_ylabel('Enthalpy (eV)')
    ax.set_xlabel('Pressure (eV/$\AA^3$)')
    plt.savefig('enthalpy_vs_pressure.png')
    return None, {}
    