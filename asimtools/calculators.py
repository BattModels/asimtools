'''
Tools for loading and returning ASE calculator objects for use in simulations
'''
import importlib
import logging
from typing import Dict, Optional
from asimtools.utils import get_calc_input

# pylint: disable=import-outside-toplevel
# pylint: disable=import-error

def load_calc(
    calc_id: str = None,
    calc_input: Optional[Dict] = None
):
    """Loads a calculator using given calc_id or calc_input

    :param calc_id: ID/key to use to load calculator from the supplied/global\
        calc_input file, defaults to None
    :type calc_id: str, optional
    :param calc_input: calc_input dictionary for a single calculator\
        , defaults to None
    :type calc_input: Optional[Dict], optional
    :return: ASE calculator instance
    :rtype: :class:`ase.calculators.calculators.Calculator`
    """
    if calc_id is not None:
        if calc_input is None:
            calc_input = get_calc_input()
        try:
            calc_params = calc_input[calc_id]
        except KeyError as exc:
            msg = f'Calculator with calc_id: {calc_id} not found in'
            msg += f'calc_input {calc_input}'
            raise KeyError(msg) from exc
        except AttributeError as exc:
            raise AttributeError('No calc_input found') from exc
    name = calc_params.get('name', None)

    if 'module' in calc_params:
        loader = load_ase_calc
    else:
        try:
            loader = external_calcs[name]
        except KeyError as exc:
            msg = 'Provide ASE module for calc or use'
            msg += f' one of {external_calcs.keys()}'
            raise KeyError(msg) from exc

    calc = loader(calc_params=calc_params)
    label = calc_params.get('label', name)
    calc.label = label
    logging.debug('Loaded calculator %s', calc_id)
    return calc


def load_nequip(calc_params):
    """Load NequIP or Allegro calculator

    :param calc_params: args to pass to loader
    :type calc_params: Dict
    :return: NequIP ase calculator
    :rtype: :class:`nequip.ase.nequip_calculator.NequIPCalculator`
    """
    from nequip.ase.nequip_calculator import NequIPCalculator
    if calc_params.get('float64', False):
        import torch
        torch.set_default_dtype(torch.float64)
    try:
        calc = NequIPCalculator.from_deployed_model(**calc_params['args'])
    except Exception as exc:
        msg = f"Failed to load NequIP with parameters:\n {calc_params}"
        raise RuntimeError(msg) from exc

    return calc


def load_dp(calc_params):
    """Load Deep Potential Calculator

    :param calc_params: args to pass to loader
    :type calc_params: Dict
    :return: DP calculator
    :rtype: :class:`deepmd.calculator.DP`
    """
    from deepmd.calculator import DP

    try:
        calc = DP(**calc_params['args'])
    except Exception:
        print(f"Failed to load DP with parameters:\n {calc_params}")
        raise

    return calc


def load_ase_calc(calc_params):
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
    'NequIP': load_nequip,
    'Allegro': load_nequip,
    'DeepPotential': load_dp,
}
