'''
Tools for loading and returning ASE calculator objects for use in simulations
'''

#pylint: disable=import-outside-toplevel
#pylint: disable=import-error

# def load_generic(calc_params):
#     ''' Loads a generic ASE calculator which follows the standard format '''
#     from ase.calculators.eam import EAM

def load_calc(calc_params):
    ''' Finds the correct loader and load the calc '''
    calcs = [calcn for calcn in calc_dict]
    name = calc_params.get('name', False)
    assert name, f'{name} calculator is not implemented. Calcs are {calcs}'
    loader = calc_dict.get(name)

    calc_input = calc_params.get('calc_input', {})
    # assert calc_input, 'Provide calc_input specific to chosen calculator'

    calc = loader(calc_input=calc_input)
    label = calc_params.get('label', name)
    calc.label = label

    return calc


def _load_nequip(name, calc_input):
    ''' Loads a NequIP or Allegro Calculator '''
    from nequip.ase.nequip_calculator import NequIPCalculator

    try:
        calc = NequIPCalculator.from_deployed_model(**calc_input)
    except Exception:
        print(f"Failed to load {name} with parameters:\n {calc_input}")
        raise

    return calc


def _load_quantumespresso(calc_input):
    ''' Loads QE calculator '''
    raise NotImplementedError


def load_lj(calc_input):
    ''' Loads a Lennard Jones Calculator '''
    from ase.calculators.lj import LennardJones

    try:
        calc = LennardJones(**calc_input)
    except Exception:
        print(f"Failed to load calculator with parameters:\n {calc_input}")
        raise

    return calc

calc_dict = {
    # 'NequIP': load_nequip,
    'LennardJones': load_lj,
}
