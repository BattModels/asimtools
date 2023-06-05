'''
Tools for loading and returning ASE calculator objects for use in simulations
'''

import importlib
# pylint: disable=import-outside-toplevel
# pylint: disable=import-error

# def load_generic(calc_params):
#     ''' Loads a generic ASE calculator which follows the standard format '''
#     from ase.calculators.eam import EAM

def load_calc(calc_params):
    ''' Finds the correct loader and load the calc '''

    name = calc_params.get('name', None)
    if 'module' in calc_params:
        loader = _load_ase_calc
    else:
        try:
            loader = external_calcs[name]
        except KeyError:
            print(f'Provide ASE module for calc or \
                one of {external_calcs.keys()}')
            raise

    calc = loader(calc_params=calc_params)
    label = calc_params.get('label', name)
    calc.label = label

    return calc


def _load_nequip(calc_params):
    ''' Loads a NequIP or Allegro Calculator '''
    from nequip.ase.nequip_calculator import NequIPCalculator

    try:
        calc = NequIPCalculator.from_deployed_model(**calc_params)
    except Exception:
        print(f"Failed to load NequIP with parameters:\n {calc_params}")
        raise

    return calc

def _load_dp(calc_params):
    ''' Loads a DP Calculator '''
    from deepmd.calculator import DP

    try:
        calc = DP(**calc_params)
    except Exception:
        print(f"Failed to load DP with parameters:\n {calc_params}")
        raise

    return calc


def _load_quantumespresso(calc_params):
    ''' Loads QE calculator and catches common errors'''
    from ase.calculators.espresso import Espresso

    # Perform tests on input
    necessary_params = ['input_data', 'pseudopotentials', 'command']
    print(calc_params)
    for param in necessary_params:
        assert param in calc_params, f'{param} parameter missing for calc'

    input_data = calc_params.get('input_data', False)
    assert input_data, 'Provide input_data for QE'
    necessary_input_data = ['ecutwfc']
    for data in necessary_input_data:
        assert data in calc_params['input_data'], \
            f'{param} parameter missing for QE'

    input_data['calculation'] = input_data.get('calculation', 100)
    if input_data['calculation'] != 'scf':
        cell_dynamics = input_data.get('cell_dynamics')
        ion_dynamics = input_data.get('ion_dynamics', False)
        assert cell_dynamics or ion_dynamics, \
            f'{input_data["calculation"]} calculation might need dynamics'

    kspacing = calc_params.get('kspacing', None)
    kpts = calc_params.get('kpts', None)
    assert (not (kpts is None and kspacing is None)) or \
        (not ((kpts is not None) and (kspacing is not None))), \
        'Provide exactly one of kpts or kspacing'

    calc = Espresso(
        pseudopotentials=calc_params['pseudopotentials'],
        tstress=True,
        tprnfor=True,
        kpts=kpts,
        kspacing=kspacing,
        input_data=calc_params['input_data'],
        command=calc_params['command']
    )

    return calc


def _load_lj(calc_params):
    ''' Loads a Lennard Jones Calculator '''
    from ase.calculators.lj import LennardJones

    try:
        calc = LennardJones(**calc_params)
    except Exception:
        print(f"Failed to load calculator with parameters:\n {calc_params}")
        raise

    return calc

def _load_ase_calc(calc_params):
    ''' Load any builtin ASE calculator '''
    module_name = calc_params.get('module', '')
    try:
        module = importlib.import_module(module_name)
    except:
        print("Check calc module")
        raise
    name = calc_params.get('name', None)
    try:
        calc_class = getattr(module, name)
    except:
        print("Check calc name")
        raise
    calc_args = calc_params.get('args', {})
    try:
        calc = calc_class(**calc_args)
    except:
        print("Check calc args")
        raise
    return calc

external_calcs = {
    'NequIP': _load_nequip,
    'Allegro': _load_nequip,
    # 'QuantumEspresso': _load_quantumespresso,
    'DeepPotential': _load_dp,
}
