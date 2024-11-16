# Change Log
All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).
 
## [develop] - 2024-11-13
 
### Added
- ASIMMODULE: benchmarking/distribution
- ASIMMODULE: data/collect_images asimmodule
- CALCULATOR: OMAT24 calculator
- VASP examples
- Jobs submitted as an array can now be grouped together in slurm jobs. See
  distributed examples
 
### Changed
- We now append rather than rewrite to stderr.txt and stdout.txt 
- Now forcing natural sorting of files in workflows file-based arrays were
  confusing

### Fixed
- Bugs in phonopy asimmodule that would give incorrectly sorted EOS energies
- If an atoms object had unknown keys, writing the output extxyz file would
  throw an error, now all the ASE sanctioned properties are written
