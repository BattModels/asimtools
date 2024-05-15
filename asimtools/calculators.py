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

    https://github.com/mir-group/nequip/tree/main

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
        calc = NequIPCalculator.from_deployed_model(**calc_params.get('args', {}))
    except Exception as exc:
        msg = f"Failed to load NequIP with parameters:\n {calc_params}"
        raise RuntimeError(msg) from exc

    return calc


def load_dp(calc_params):
    """Load Deep Potential Calculator

    https://docs.deepmodeling.com/projects/deepmd/en/master/

    :param calc_params: args to pass to loader
    :type calc_params: Dict
    :return: DP calculator
    :rtype: :class:`deepmd.calculator.DP`
    """
    from deepmd.calculator import DP

    try:
        calc = DP(**calc_params['args'])
    except Exception:
        logging.error("Failed to load DP with parameters:\n %s", calc_params)
        raise

    return calc


def load_ase_calc(calc_params):
    ''' Load any builtin ASE calculator '''
    module_name = calc_params.get('module', '')
    try:
        module = importlib.import_module(module_name)
    except:
        logging.error("Check calc module and stderr")
        raise
    name = calc_params.get('name', None)
    try:
        calc_class = getattr(module, name)
    except:
        logging.error("Check calc name and stderr")
        raise
    calc_args = calc_params.get('args', {})
    try:
        calc = calc_class(**calc_args)
    except:
        logging.error("Check calc args and stderr")
        raise
    return calc

def load_chgnet(calc_params):
    """Load CHGNet Calculator

    https://chgnet.lbl.gov/#tutorials-and-docs

    :param calc_params: args to pass to loader
    :type calc_params: Dict
    :return: CHGNet calculator
    :rtype: :class:`chgnet.model.dynamics.CHGNetCalculator`
    """
    from chgnet.model.dynamics import CHGNetCalculator
    if calc_params['args'].get('from_file', False):
        calc_params['args'].pop('from_file')
        try:
            calc = CHGNetCalculator.from_file(**calc_params['args'])
        except Exception:
            logging.error("Failed to load CHGNet from file with parameters:\
                \n %s", calc_params)
            raise
    else:
        try:
            calc = CHGNetCalculator(**calc_params['args'])
        except Exception:
            logging.error("Failed to load CHGNet with parameters:\n %s", \
                calc_params)
            raise

    return calc

def load_mace(calc_params):
    """Load MACE Calculator

    https://github.com/ACEsuit/mace?tab=readme-ov-file

    :param calc_params: args to pass to loader
    :type calc_params: Dict
    :return: MACE calculator
    :rtype: :class:`mace.calculators.mace_mp`
    """
    from mace.calculators import mace_mp

    try:
        calc = mace_mp(**calc_params['args'])
    except Exception:
        logging.error("Failed to load MACE with parameters:\n %s", calc_params)
        raise

    return calc

def load_mace_off(calc_params):
    """Load MACE Calculator

    https://github.com/ACEsuit/mace?tab=readme-ov-file

    :param calc_params: args to pass to loader
    :type calc_params: Dict
    :return: MACE-OFF calculator
    :rtype: :class:`mace.calculators.mace_mp`
    """
    from mace.calculators import mace_off

    try:
        calc = mace_off(**calc_params['args'])
    except Exception:
        logging.error("Failed to load MACE-OFF with parameters:\n %s", calc_params)
        raise

    return calc

def load_espresso_profile(calc_params):
    """Load Qunatum Espresso Calculator for ASE>3.22.1. If using older versions
    such as the ones available on PyPI or conda-forge, just load it as an ASE
    calculator. Until the new ASE version becomes an official release, we will
    have to have both for compatibility. The interface however remains that of 
    ASE <=3.22.1 within ASIMTools for consistency using the `command` keyword

    https://wiki.fysik.dtu.dk/ase/releasenotes.html

    :param calc_params: args to pass to loader
    :type calc_params: Dict
    :return: Espresso calculator
    :rtype: :class:`ase.calculators.espresso.Espresso`
    """
    from ase.calculators.espresso import Espresso, EspressoProfile

    if 'command' in calc_params['args']:
        command = calc_params['args'].pop('command')
        command = command.split()
        progind = command.index('pw.x')
        argv = command[:progind+1]
    else:
        argv = ['pw.x']

    try:
        calc = Espresso(
            **calc_params['args'],
            profile=EspressoProfile(argv=argv)
        )
    except Exception:
        logging.error("Failed to load MACE-OFF with parameters:\n %s", calc_params)
        raise

    return calc


def load_m3gnet(calc_params):
    """Load and M3GNet calculator

    :param calc_params: parameters to be passed to matgl.ext.ase.M3GNetCalculator. Must include a key "model" that points to the model used to instantiate the potential
    :type calc_params: Dict
    :return: M3GNet calculator
    :rtype: :class:`matgl.ext.ase.M3GNetCalculator`
    """
    from matgl.ext.ase import M3GNetCalculator
    import matgl

    model = calc_params.pop("model")
    try:
        pot = matgl.load_model(model)
        calc = M3GNetCalculator(
            pot,
            **calc_params,
        )
    except Exception:
        logging.error("Failed to load M3GNet with parameters:\n %s", calc_params)
        raise

    return calc

external_calcs = {
    'NequIP': load_nequip,
    'Allegro': load_nequip,
    'DeepPotential': load_dp,
    'CHGNet': load_chgnet,
    'MACE': load_mace,
    'EspressoProfile': load_espresso_profile,
    'M3GNet': load_m3gnet,
}
