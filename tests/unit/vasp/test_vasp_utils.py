'''
Tests for asimmodules/vasp/utils.py
'''
import pytest
from unittest.mock import patch, MagicMock
from asimtools.asimmodules.vasp.utils import write_vasp_inputs

IMAGE = {'name': 'Na', 'builder': 'bulk'}


def test_write_vasp_inputs_invalid_mpset():
    ''' Unknown mpset name raises ImportError '''
    with pytest.raises(ImportError, match='Unknown mpset'):
        write_vasp_inputs(IMAGE, mpset='NotARealSet')


def test_write_vasp_inputs_with_mpset():
    ''' Valid mpset instantiates the set class and calls write_input '''
    mock_vasp_input = MagicMock()
    mock_set_cls = MagicMock(return_value=mock_vasp_input)

    with patch('pymatgen.io.vasp.sets.MPRelaxSet', mock_set_cls):
        result = write_vasp_inputs(IMAGE, mpset='MPRelaxSet')

    mock_set_cls.assert_called_once()
    mock_vasp_input.write_input.assert_called_once_with("./")
    assert result is mock_vasp_input


def test_write_vasp_inputs_with_mpset_forwards_incar_settings():
    ''' user_incar_settings is forwarded to the set constructor '''
    mock_vasp_input = MagicMock()
    mock_set_cls = MagicMock(return_value=mock_vasp_input)
    incar_settings = {'ENCUT': 520, 'KSPACING': 0.22}

    with patch('pymatgen.io.vasp.sets.MPRelaxSet', mock_set_cls):
        write_vasp_inputs(
            IMAGE,
            mpset='MPRelaxSet',
            user_incar_settings=incar_settings,
        )

    _, kwargs = mock_set_cls.call_args
    assert kwargs['user_incar_settings'] == incar_settings


def test_write_vasp_inputs_with_prev_calc():
    ''' prev_calc triggers from_prev_calc instead of direct construction '''
    mock_vasp_input = MagicMock()
    mock_set_cls = MagicMock()
    mock_set_cls.from_prev_calc.return_value = mock_vasp_input

    with patch('pymatgen.io.vasp.sets.MPRelaxSet', mock_set_cls):
        result = write_vasp_inputs(
            IMAGE,
            mpset='MPRelaxSet',
            prev_calc='/path/to/prev',
        )

    mock_set_cls.from_prev_calc.assert_called_once()
    mock_set_cls.assert_not_called()
    mock_vasp_input.write_input.assert_called_once_with("./")
    assert result is mock_vasp_input


def test_write_vasp_inputs_no_mpset():
    ''' No mpset builds VaspInput manually and calls check_params on Incar '''
    mock_incar = MagicMock()
    mock_vasp_input = MagicMock()

    with patch('asimtools.asimmodules.vasp.utils.Incar', return_value=mock_incar) as mock_incar_cls, \
         patch('asimtools.asimmodules.vasp.utils.Poscar'), \
         patch('asimtools.asimmodules.vasp.utils.VaspInput', return_value=mock_vasp_input):

        result = write_vasp_inputs(IMAGE, user_incar_settings={'ENCUT': 400})

    mock_incar_cls.assert_called_once_with({'ENCUT': 400})
    mock_incar.check_params.assert_called_once()
    mock_vasp_input.write_input.assert_called_once_with("./")
    assert result is mock_vasp_input


def test_write_vasp_inputs_no_mpset_with_kpoints():
    ''' user_kpoints_settings causes a Kpoints object to be constructed '''
    mock_kpoints = MagicMock()

    with patch('asimtools.asimmodules.vasp.utils.Incar', return_value=MagicMock()), \
         patch('asimtools.asimmodules.vasp.utils.Poscar'), \
         patch('asimtools.asimmodules.vasp.utils.Kpoints', return_value=mock_kpoints) as mock_kpoints_cls, \
         patch('asimtools.asimmodules.vasp.utils.VaspInput', return_value=MagicMock()) as mock_vasp_input_cls:

        write_vasp_inputs(IMAGE, user_kpoints_settings={'kpts': [4, 4, 4]})

    mock_kpoints_cls.assert_called_once_with({'kpts': [4, 4, 4]})
    _, kwargs = mock_vasp_input_cls.call_args
    assert kwargs['kpoints'] is mock_kpoints


def test_write_vasp_inputs_no_mpset_no_kpoints():
    ''' No user_kpoints_settings passes kpoints=None to VaspInput '''
    with patch('asimtools.asimmodules.vasp.utils.Incar', return_value=MagicMock()), \
         patch('asimtools.asimmodules.vasp.utils.Poscar'), \
         patch('asimtools.asimmodules.vasp.utils.VaspInput', return_value=MagicMock()) as mock_vasp_input_cls:

        write_vasp_inputs(IMAGE)

    _, kwargs = mock_vasp_input_cls.call_args
    assert kwargs['kpoints'] is None


def test_write_vasp_inputs_vaspinput_kwargs_forwarded():
    ''' vaspinput_kwargs are unpacked into the set constructor '''
    mock_vasp_input = MagicMock()
    mock_set_cls = MagicMock(return_value=mock_vasp_input)
    extra = {'optional_structure_charge': 0}

    with patch('pymatgen.io.vasp.sets.MPRelaxSet', mock_set_cls):
        write_vasp_inputs(IMAGE, mpset='MPRelaxSet', vaspinput_kwargs=extra)

    _, kwargs = mock_set_cls.call_args
    assert kwargs.get('optional_structure_charge') == 0
