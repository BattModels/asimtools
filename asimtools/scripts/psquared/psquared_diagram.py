'''.Xauthority'''

from typing import Dict, Tuple
from glob import glob
import numpy as np
from scipy import interpolate
from ase.units import GPa
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry, Element, PDPlotter
from pymatgen.core import Composition
from asimtools.utils import write_csv_from_dict
matplotlib.use('Agg')

def plot_phd(x, y, Z, qhull_formulas, fname='phase_diagram.png', log=True):
    X, Y = np.meshgrid(x,y)
    fig, ax = plt.subplots()
    cmap = plt.cm.get_cmap('Paired', len(qhull_formulas)) 
    c = ax.pcolor(X, Y, Z, cmap=cmap)

    plt.legend(
        [mpatches.Patch(color=cmap(f)) for f in range(len(qhull_formulas))],
        [f'{qhull_formula}' for qhull_formula in qhull_formulas]
    )
    if log:
        plt.xscale('log')
    plt.savefig(fname)

# @branch
def psquared_diagram(
    gp_files_pattern: str,
    pH: float,
    H_ref_csv: str,
    Temp: float = 300,
    plogspace: Tuple[float,float,float] = (-5,-2,10),
    ulim: Tuple[float,float] = (-1.1,0.2),
    nUs: int = 200,
) -> Tuple[None,Dict]:
    '''
    Generating P2 phase diagrams given the Gibbs energy versus pressure. Note
    that this works with either formation energies or without since if you give
    formation energies the references will be zero

    This assumes that the Gibbs free energies provided are the minimum gibbs
    for all polymorphs of that composition
    '''
    
    # For comp in comps directory:
    gp_csvs = glob(gp_files_pattern)
    gp_func_dict = {}
    for gp_csv in gp_csvs:
        with open(gp_csv, 'r', encoding='utf-8') as f:
            formula = f.readline().split(',')[0].split(':')[-1]
            # To match pymatgen format
            formula = Composition(formula).formula

        gp_data = np.genfromtxt(gp_csv, delimiter=',')
        pvals = gp_data[:,0]
        gvals = gp_data[:,1]
        print(formula, 'range:', np.min(gp_data[:,0]), np.max(gp_data[:,0]))
        # Make interpolation function for each comp
        gp_func = interpolate.interp1d(pvals, gvals, kind='cubic')
        gp_func_dict[formula] = gp_func

    # Get H reference data, assuming it is given per atom here
    href_data = np.genfromtxt(H_ref_csv, delimiter=',')
    gp_func = interpolate.interp1d(
        href_data[:,0]*GPa, href_data[:,1], kind='cubic'
    )
    print('hprange:', np.min(href_data[:,0]*GPa), np.max(href_data[:,0]*GPa))
    gp_func_dict['H2'] = gp_func

    # # Get the key for the pure element that isn't hydrogen
    # elems = [key for key in gp_func_dict if 'H' not in key]
    # assert len(elems) == 1, 'Provide only 1 pure element reference'
    # M_key = elems[0]
    # # Number of M atoms in data for ref
    # refnM = Composition(M_key).num_atoms

    ps = np.logspace(*plogspace) * GPa
    formulas = list(gp_func_dict.keys())

    p2_grid = np.zeros([len(ps), nUs])
    for pp, pressure in enumerate(ps):
        # GM = gp_func_dict[M_key]
        # GH = gp_func_dict['Hs']
        G_entries = []
        for formula in formulas:
            print(formula)
            try:
                GMHn = gp_func_dict[formula](pressure)
            except:
                print('Error:', formula, pressure)
                raise
            entry = PDEntry(formula, GMHn)
            G_entries.append(entry)
        pd_formulas = []
        Gf_energies = []
        phase_diagram = PhaseDiagram(G_entries)
        for entry in G_entries:
            pd_formulas.append(entry.composition.formula)
            Gf = phase_diagram.get_form_energy(entry)
            Gf_energies.append(Gf)

            # Try using Gf/formula unit
            # formula, factor = entry.composition.get_integer_formula_and_factor()
            # pd_formulas.append(formula)
            # Gf_energies.append(phase_diagram.get_form_energy(entry) / factor)

        # For now use Pinwen's method of setting GH=1e5
        H_ind = pd_formulas.index('H2')
        pd_formulas.pop(H_ind)
        Gf_energies.pop(H_ind)
        # import pdb; pdb.set_trace()

        print(f'p: {pressure/GPa}\nformulas = {formulas}\nformulas = {pd_formulas}\nenergies = {Gf_energies}')
        # pb_grid, grid_formulas = convhull_pourbaix(
        #     {},
        #     formulas=pd_formulas[:] + ['H2'],
        #     energies=Gf_energies[:] + [1e5],
        #     nUs=nUs,
        #     npHs=2,
        #     ulim=ulim,
        #     # pHlim=[pH,pH+0.001], #No need for full diagram, just a slice
        #     pHlim=[-2,14],
        #     T=Temp,
        #     plot=True,
        #     # plot=False,
        #     plot_fname=f'pb/p{pressure:.2E}_pourbaix.png'
        # )
        pb_grid, grid_formulas = Pourbaix_H_2(
            formula=pd_formulas[:],
            E=Gf_energies[:],
            nu=nUs,
            nph=300,
            ulim=(-1.1,0.2),
            phlim=(-2,14),
            plot=False,
            plot_fname=f'pb/p{pressure:.2E}_pourbaix.png',
        )
        # pb_slice = pb_grid[:,0] # Get only the slice for the pH we want
        # print('shape', pb_grid)
        # print(grid_formulas)
        # p2_grid[pp] = pb_slice

        ph_ind = int((pH-(-2))/16 - 2)
        p2_grid[pp] = pb_grid[:-1,ph_ind]
        

    p2_grid = np.transpose(p2_grid)
    np.savetxt('p2_diagram.csv', p2_grid, delimiter=',')
    write_csv_from_dict('formulas.txt', {'formulas': formulas})
    plot_phd(ps / GPa, np.linspace(*ulim, nUs), p2_grid, grid_formulas, fname='p2_diagram.png')
    return {}

def convhull_pourbaix(
    formulas,
    energies,
    nUs=200,
    npHs=300,
    ulim=(-1.1,0.2),
    pHlim=(-2,14),
    T=300,
    plot=False,
    plot_fname='pourbaix_diagram.png',
    **kwargs
):
    ''' 
    Generate Pourbaix using structures with intercept method implemented
    by Pinwen
    '''

    assert len(energies) == len(formulas), 'energies and formulae'
    assert 'H' not in formulas[0], '1st one should be pure without H'
    # assert len(set(formulas[-1])) == 1, f'{set(formulas[-1])} Last one should be H'

    Us = np.linspace(*ulim, nUs)
    pHs = np.linspace(*pHlim, npHs)

    entries = [PDEntry(formulas[i], energies[i]) for i in range(len(formulas))]
    phd =  PhaseDiagram(entries)
    # plot = PDPlotter(phd, backend='matplotlib')
    # plot.get_plot()
    # plt.savefig(f'convhulls/hull{np.random.randint(1000)}.png')
    # Get compositions on the hull and their hydrogen fractions
    stable_entries = [
        entry for entry in phd.qhull_entries if \
        phd.get_decomp_and_e_above_hull(entry,allow_negative=True)[1] == 0
    ][:-1]
    qhull_es_per_atom = np.array([entry.energy_per_atom for entry in stable_entries])
    qhull_formulas = np.array([entry.composition.reduced_formula for entry in stable_entries])
    # print(qhull_formulas)
    # print(qhull_es_per_atom)
    xs = np.array([entry.composition.get_atomic_fraction(Element('H')) for entry in stable_entries])

    # Define elemental references and some useful constants
    # M_energy = energies[0]
    # H_energy = energies[-1]
    kB = 8.617333e-5 #eV/K
    gamma = kB * T * np.log(10) # This is just a constant we simplify here

    # Sort compositions by H fraction
    sort_ind = xs.argsort()
    xs = xs[sort_ind]
    gfMHns = qhull_es_per_atom[sort_ind]
    qhull_formulas = qhull_formulas[sort_ind]

    # Get the intercept at zero U and lnpH=0
    muMHn_intercepts = gfMHns[1:] + np.multiply(
        np.divide(gfMHns[1:]-gfMHns[:-1], (xs[1:]-xs[:-1])), (1-xs[1:]),
    )

    # Generate the phase grid by setting the phase with the highest
    # hydrogen content as the most stable once the chemical potential
    # of H has changed enough such that DeltaG for hydride is less than 0
    phase_grid = -1 * np.ones([nUs, npHs])
    for uu, U in enumerate(Us):
        for MHn_index, muMHn in enumerate(muMHn_intercepts):
            max_stable_pH = -(U + muMHn) / gamma

            for pp, pH in enumerate(pHs):
                if pH <= max_stable_pH:
                    phase_grid[uu,pp] = MHn_index

    if plot:
        plot_phd(pHs, Us, phase_grid, qhull_formulas, fname=plot_fname, log=False)
    return phase_grid, qhull_formulas

def Pourbaix_H_2(
    formula,
    E,
    nu=200,
    nph=300,
    ulim=(-1.1,0.2),
    phlim=(-2,14),
    plot=False,
    plot_fname='pourbaix.png',
    **kwargs
):
    # Unnecessary reloading of module
    from pymatgen.analysis.phase_diagram import PhaseDiagram
    
    # Unnecessary relabelling of variable
    comp=formula
    
    # Unnecessary relabeling of variable
    y=E

    # GEt numuber of atoms for each entry
    nat=[Composition(i).num_atoms for i in comp]

    # Individual entry to PD, Set H energy high so that all hydrides show up in phase diagram
    entry=[PDEntry(comp[i],y[i]) for i in range(len(y))]+[PDEntry('H',1e5)]
    
    # Energy above convex hull for each entry
    stability=[PhaseDiagram(entry).get_decomp_and_e_above_hull(i,allow_negative=True)[1] for i in entry[:-1]]

    phd =  PhaseDiagram(entry)
    # plot = PDPlotter(phd, backend='matplotlib')
    # plot.get_plot()
    # plt.savefig(f'convhulls/hull{np.random.randint(1000)}.png')
    
    # Compositions on the convex hull
    comp_convex=[comp[i] for i in range(len(comp)) if stability[i]==0]

    # energies per atom on the convex hull. Shouldn't this be per formula unit?
    y_convex=[y[i]/nat[i] for i in range(len(y)) if stability[i]==0]

    # Fraction of hydrogen
    x=[1-float(i.composition.fractional_composition.formula.split(' ')[0].strip(comp[0])) for i in entry[:-1]]

    # Hydrogen fractions on convex hull
    x_convex=[x[i] for i in range(len(x)) if stability[i]==0]

    # What is mutr? mutr is the chemical potential change below which(as it gets more negative) the compound at i+1 is no longer on the convex hull
    mutr=['' for i in range(len(comp_convex)-1)]
    for i in range(len(comp_convex)-1):
        mutr[i]=y_convex[i]+(1-x_convex[i])*(y_convex[i]-y_convex[i+1])/(x_convex[i]-x_convex[i+1])

    # print(y_convex, mutr)
    # T=300 K, U=-0.0596*pH-mutr
    phmin=phlim[0];phmax=phlim[1];umin=ulim[0];umax=ulim[1]
    a=-np.ones(shape=(nu+1,nph+1))
    for i in range(nu+1):
        u=umin+i/nu*(umax-umin) # Equivalent to np.linspace(umin, umax, nu)
        for j in range(len(mutr)+1):
            # This essentially finds the ranges of pH for which the particular compound
            # is stable
            if j==0:
                # This is the metal end of the hull therefore we have maximum alkalinity
                ph1=phmax
            else:
                # Compare the chemical potential minus applied electric potential, if 
                ph1=min((-u-mutr[j-1])/0.0596,phmax)
            # print(j, ph1, mutr, u)
            if j==len(mutr):
                # This is the hydrogen end of the hull therefore we have max acidity
                ph2=phmin
            else:
                ph2=max((-u-mutr[j])/0.0596,phmin)
            for k in range(
                int(np.ceil((ph2-phmin)/(phmax-phmin)*nph)),
                int(np.floor((ph1-phmin)/(phmax-phmin)*nph))+1
            ):
                a[i,k]=j
    if plot:
        phs = np.linspace(*phlim, nph+1)
        us = np.linspace(*ulim, nu+1)
        plot_phd(phs, us, a, comp_convex, fname=plot_fname, log=False)
    return a,comp_convex
