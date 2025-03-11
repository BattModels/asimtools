# Change Log
All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [develop] - 2025-2-14
 
### Added
- Can now use placehodlers in sim_array where the array_values replace part of
the string at the given key_sequence
- Can now specify whether to write velocities in lammps
- Calculator: Now added DFTD3 calculator using the ASE interface, works with any calculator. Two things of note: 1. You need to install DFTD3 from 
https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3/get_dft-d3. 2. Some
calculators which return a 3x3 matrix for stress will break. One can modify
ASE source for this as ASIMTools can't go into the calculator code.
 
### Changed
- VASP interface changed to align more with pymatgen
- asimtools.utils.write_yaml now stops sorting keys to help with readability of 
written yamls

### Fixed
- Minor bugs in geometry optimizations
- Updated EspressoProfile calculator to match ASE 3.25.0b1

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
