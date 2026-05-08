'''
Tests for calculators.py
'''
import pytest
from ase.calculators.emt import EMT
from asimtools.calculators import load_calc, load_ase_calc

EMT_PARAMS = {
    'name': 'EMT',
    'module': 'ase.calculators.emt',
    'args': {},
}

CALC_INPUT = {'emt': EMT_PARAMS}


# ── load_ase_calc ─────────────────────────────────────────────────────────────

def test_load_ase_calc_returns_correct_type():
    ''' load_ase_calc with EMT params returns an EMT instance '''
    calc = load_ase_calc(EMT_PARAMS)
    assert isinstance(calc, EMT)


def test_load_ase_calc_with_args():
    ''' load_ase_calc passes args to the calculator constructor '''
    params = {
        'name': 'LennardJones',
        'module': 'ase.calculators.lj',
        'args': {'epsilon': 2.0, 'sigma': 1.5},
    }
    from ase.calculators.lj import LennardJones
    calc = load_ase_calc(params)
    assert isinstance(calc, LennardJones)
    assert calc.parameters['epsilon'] == 2.0
    assert calc.parameters['sigma'] == 1.5


def test_load_ase_calc_bad_module():
    ''' load_ase_calc raises when module cannot be imported '''
    params = {'name': 'EMT', 'module': 'nonexistent.module', 'args': {}}
    with pytest.raises(ModuleNotFoundError):
        load_ase_calc(params)


def test_load_ase_calc_bad_name():
    ''' load_ase_calc raises when class name is not found in module '''
    params = {'name': 'NotAClass', 'module': 'ase.calculators.emt', 'args': {}}
    with pytest.raises(AttributeError):
        load_ase_calc(params)


# ── load_calc: new calculator dict interface ──────────────────────────────────

def test_load_calc_calculator_with_calc_params():
    ''' calculator={"calc_params": ...} loads without a calc_input lookup '''
    calc = load_calc(calculator={'calc_params': EMT_PARAMS})
    assert isinstance(calc, EMT)


def test_load_calc_calculator_with_calc_id():
    ''' calculator={"calc_id": ...} looks up params from supplied calc_input '''
    calc = load_calc(
        calculator={'calc_id': 'emt'},
        calc_input=CALC_INPUT,
    )
    assert isinstance(calc, EMT)


def test_load_calc_calculator_with_calc_id_global(tmp_path, monkeypatch):
    ''' calculator={"calc_id": ...} falls back to global calc_input env var '''
    from asimtools.utils import write_yaml
    calc_input_file = tmp_path / 'calc_input.yaml'
    write_yaml(calc_input_file, CALC_INPUT)
    monkeypatch.setenv('ASIMTOOLS_CALC_INPUT', str(calc_input_file))
    calc = load_calc(calculator={'calc_id': 'emt'})
    assert isinstance(calc, EMT)


def test_load_calc_calculator_sets_label():
    ''' load_calc sets calc.label from calc_params name when no label given '''
    calc = load_calc(calculator={'calc_params': EMT_PARAMS})
    assert calc.label == 'EMT'


def test_load_calc_calculator_sets_custom_label():
    ''' load_calc sets calc.label from explicit label key in calc_params '''
    params = dict(EMT_PARAMS, label='my_emt')
    calc = load_calc(calculator={'calc_params': params})
    assert calc.label == 'my_emt'


# ── load_calc: legacy interface (backward compatibility) ──────────────────────

def test_load_calc_legacy_calc_params():
    ''' Legacy calc_params kwarg still works '''
    calc = load_calc(calc_params=EMT_PARAMS)
    assert isinstance(calc, EMT)


def test_load_calc_legacy_calc_id():
    ''' Legacy calc_id kwarg still works with explicit calc_input '''
    calc = load_calc(calc_id='emt', calc_input=CALC_INPUT)
    assert isinstance(calc, EMT)


# ── load_calc: error cases ────────────────────────────────────────────────────

def test_load_calc_no_args_raises():
    ''' load_calc raises AssertionError when called with no identifying args '''
    with pytest.raises(AssertionError):
        load_calc()


def test_load_calc_missing_calc_id_raises():
    ''' load_calc raises KeyError when calc_id is not in calc_input '''
    with pytest.raises(KeyError):
        load_calc(calc_id='missing', calc_input=CALC_INPUT)


def test_load_calc_no_module_or_external_name_raises():
    ''' load_calc raises KeyError when calc_params has unknown name and no module '''
    params = {'name': 'UnknownCalc', 'args': {}}
    with pytest.raises(KeyError):
        load_calc(calc_params=params)


def test_load_calc_calculator_takes_precedence_over_legacy():
    ''' calculator kwarg overrides legacy calc_params when both provided '''
    other_params = {'name': 'LennardJones', 'module': 'ase.calculators.lj', 'args': {}}
    calc = load_calc(
        calculator={'calc_params': EMT_PARAMS},
        calc_params=other_params,
    )
    assert isinstance(calc, EMT)
