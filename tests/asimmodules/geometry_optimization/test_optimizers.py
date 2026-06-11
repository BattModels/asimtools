"""
Tests for the geometry_optimization asimmodules using ASE's EMT calculator
on copper.
"""
import pytest
from asimtools.utils import get_atoms
from asimtools.asimmodules.geometry_optimization.atom_relax import atom_relax
from asimtools.asimmodules.geometry_optimization.cell_relax import cell_relax
from asimtools.asimmodules.geometry_optimization.optimize import optimize
from asimtools.asimmodules.geometry_optimization.symmetric_cell_relax import (
    symmetric_cell_relax,
)
from asimtools.asimmodules.geometry_optimization.ase_cubic_eos_optimization import (
    ase_cubic_eos_optimization,
)

# EMT relaxation backends that act on the cell (and report a max stress)
CELL_RELAXERS = {
    'cell_relax': cell_relax,
    'optimize': optimize,
    'symmetric_cell_relax': symmetric_cell_relax,
}


@pytest.fixture
def emt_calculator():
    """Calculator spec using ASE's EMT, supplied directly via calc_params"""
    return {
        'calc_params': {
            'name': 'EMT',
            'module': 'ase.calculators.emt',
            'args': {},
        }
    }


@pytest.fixture
def cu_image():
    """Off-equilibrium bulk fcc copper so cell relaxations have work to do.

    EMT's equilibrium lattice constant for Cu is ~3.6 Ang, so starting at
    3.7 Ang leaves a finite stress to relax away.
    """
    return {
        'builder': 'bulk',
        'name': 'Cu',
        'crystalstructure': 'fcc',
        'a': 3.7,
        'cubic': True,
    }


@pytest.fixture
def rattled_cu_image():
    """Bulk copper with displaced atoms so atom_relax has forces to minimize"""
    return {
        'builder': 'bulk',
        'name': 'Cu',
        'crystalstructure': 'fcc',
        'a': 3.6,
        'cubic': True,
        'repeat': (2, 2, 2),
        'rattle_stdev': 0.1,
    }


def _emt_energy(image):
    """Potential energy of an image fixture using a fresh EMT calculator"""
    from ase.calculators.emt import EMT
    atoms = get_atoms(**image)
    atoms.calc = EMT()
    return atoms.get_potential_energy()


def test_atom_relax(emt_calculator, rattled_cu_image, tmp_path, monkeypatch):
    """atom_relax should converge below fmax and lower the energy"""
    monkeypatch.chdir(tmp_path)
    fmax = 0.05
    initial_energy = _emt_energy(rattled_cu_image)

    results = atom_relax(
        emt_calculator,
        rattled_cu_image,
        fmax=fmax,
        steps=200,
        prefix='test',
    )

    assert results['converged'] is True
    # final_fmax is a true per-atom force magnitude, so it must meet fmax
    assert results['final_fmax'] <= fmax
    assert results['energy'] < initial_energy
    # atom_relax does not move the cell, so no stress is reported
    assert 'final_smax' not in results

    assert (tmp_path / 'test_atom_relax.traj').exists()
    assert (tmp_path / 'test_image_output.xyz').exists()
    assert results['files'] == {
        'image': 'test_image_output.xyz',
        'traj': 'test_atom_relax.traj',
    }


@pytest.mark.parametrize('asimmodule_name', list(CELL_RELAXERS))
def test_cell_relaxers(emt_calculator, cu_image, asimmodule_name, tmp_path, monkeypatch):
    """cell_relax/optimize/symmetric_cell_relax should converge, report a
    finite max stress and lower the energy of an off-equilibrium cell"""
    monkeypatch.chdir(tmp_path)
    relax = CELL_RELAXERS[asimmodule_name]
    initial_energy = _emt_energy(cu_image)

    results = relax(
        emt_calculator,
        cu_image,
        fmax=0.05,
        steps=200,
        prefix='test',
    )

    assert results['converged'] is True
    assert results['energy'] < initial_energy
    # Cell relaxers report the maximum absolute stress component
    assert 'final_smax' in results
    assert results['final_smax'] >= 0.0
    assert results['final_fmax'] >= 0.0

    assert (tmp_path / f'test_{asimmodule_name}.traj').exists()
    assert (tmp_path / 'test_image_output.xyz').exists()


def test_cell_relax_reaches_equilibrium(emt_calculator, cu_image, tmp_path, monkeypatch):
    """A full cell+atom relax should drive Cu toward its EMT equilibrium
    lattice constant (~3.6 Ang, i.e. cubic cell volume ~46.7 Ang^3)"""
    monkeypatch.chdir(tmp_path)
    results = optimize(emt_calculator, cu_image, fmax=0.01, steps=200)
    relaxed = get_atoms(image_file=str(tmp_path / 'image_output.xyz'))
    # 4-atom cubic cell, equilibrium a ~ 3.6 Ang -> volume ~ 46.7 Ang^3
    assert relaxed.get_volume() == pytest.approx(46.7, abs=2.0)
    assert results['final_smax'] <= 0.05


def test_ase_cubic_eos_optimization(emt_calculator, cu_image, tmp_path, monkeypatch):
    """The EOS module should fit a positive bulk modulus and equilibrium
    volume and write its output files"""
    monkeypatch.chdir(tmp_path)
    results = ase_cubic_eos_optimization(
        emt_calculator,
        cu_image,
        npoints=5,
        plot=False,
        prefix='test',
    )

    assert results['v0'] > 0.0
    assert results['e0'] < 0.0  # EMT Cu cohesive energy is negative
    assert results['B'] > 0.0
    assert results['x0'] > 0.0

    assert (tmp_path / 'test_eos.csv').exists()
    assert (tmp_path / 'test_eos.traj').exists()
    assert (tmp_path / 'test_eos_image_output.xyz').exists()
