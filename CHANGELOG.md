# Change Log
All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [0.2.0] - 2026-04-06

### Breaking changes
- `change_dict_value`, `change_dict_values`, `expand_wildcards`: parameter `d` renamed to `dct`
- `get_str_btn`, `get_nth_label`: parameter `s` renamed to `string`
- `repeat_to_N` renamed to `repeat_to_n`
- All call sites in `asimmodules/workflows/` updated to use the new parameter names

### Added
- Pre-commit hooks: pylint (10.00/10 on `src/asimtools/` excluding asimmodules) and pytest run on every commit
- Tests for previously untested utilities: `strip_symbols`, `get_axis_lims`, `improve_plot`, `write_csv_from_dict`, `new_db`, `check_if_slurm_job_is_running`
- Tests for previously untested job functionality: `Job` getters/updaters, `DistributedJob` init and inline submit, `ChainedJob` init and last output, `load_job_from_directory`, `get_subjobs`, `load_job_tree`
- Sphinx docs now include all previously missing asimmodule packages: `active_learning`, `ase_md`, `data`, `mace`, `phonopy`, `vasp`
- Workflows documentation page with parameter tables and YAML examples for all 7 workflow modules: `sim_array`, `image_array`, `calc_array`, `distributed`, `chained`, `iterative`, `update_dependencies`
- `autodoc_mock_imports` in Sphinx config for optional heavy dependencies (`maml`, `mace`, `matgl`, `chgnet`, `phonopy`, `seekpath`)

### Changed
- Package moved to `src/` layout; `pyproject.toml` and pytest config updated accordingly
- Sphinx docs home page is now the project README
- `get_sim_input`, `get_calc_input`, `get_env_input` now return internal references instead of deepcopies
- `Job.__init__` still deepcopies inputs at construction time to preserve isolation from callers

### Fixed
- `asim_run.py`: precommands were executed twice (once discarding output, once capturing); now runs once
- `asim_execute.py`: `calc_input` was assigned twice before the conditional read
- `job.py` `start()` and `complete()`: two sequential `update_output` calls (2 reads + 2 writes) collapsed into one
- `asim_check.py`: each job tree node read `output.yaml` twice per print call; now reads once per node

### Performance
- `get_status()`: running-job status update no longer triggers a redundant `output.yaml` re-read
- `DistributedJob.__init__`: `use_slurm`/`use_sh` classification reduced from 3 O(n) passes to a single pass with early exit
- `load_job_from_directory`: replaced `exists()` + `read_yaml` (2 filesystem calls) with try/except around `read_yaml` (1 call)

## [develop] - 2025-2-14
 
### Added
- Can now specify an integer N to labels keyword add will use the N-th label for array workflows
- Can now use placehodlers in sim_array where the array_values replace part of
the string at the given key_sequence
- Can now specify whether to write velocities in lammps
- Calculator: Now added DFTD3 calculator using the ASE interface, works with any calculator. Two things of note: 1. You need to install DFTD3 from 
https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3/get_dft-d3. 2. Some
calculators which return a 3x3 matrix for stress will break. One can modify
ASE source for this as ASIMTools can't go into the calculator code.
- asim_check now also reports the job_ids
- get_atoms now allows addition of FixAtoms constraint
- output.yaml now includes hostname
 
### Changed
- VASP interface changed to align more with pymatgen
- asimtools.utils.write_yaml now stops sorting keys to help with readability of 
written yamls
- write_atoms now more universally used and recommended in asimmodules
- VASP calculation now fails if the max number of iterations is reached for ibrion=1,2,3

### Fixed
- Minor bugs in geometry optimizations
- Updated EspressoProfile calculator to match ASE 3.25.0b1
- FixSymmetry now  imported from ase.constraints

## [0.1.0] - 2024-12-27
 
### Added
- ASIMMODULE: vasp/* (For running vasp without ASE calculator and using MP sets)
- ASIMMODULE: active_learning/* (Only basic components of the workflow 
  and one step example with MACE. expect major improvements in next release)
- ASIMMODULE: benchmarking/distribution
- ASIMMODULE: data/collect_images asimmodule
- CALCULATOR: OMAT24 calculator
- VASP examples
- Jobs submitted as an array can now be grouped together in slurm jobs. See
  distributed examples
 
### Changed
- We now append rather than rewrite to stderr.txt and stdout.txt 
- Now forcing natural sorting of files in workflows file-based arrays which were
  confusing

### Fixed
- Bugs in phonopy asimmodule that would give incorrectly sorted EOS energies
- If an atoms object had unknown keys, writing the output extxyz file would
  throw an error, now only the ASE sanctioned properties are written
