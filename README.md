<!-- <img src="../figures/logo.png" alt="drawing" width="150"/> -->
<!-- [![code coverage](https://img.shields.io/codecov/c/gh/materialsproject/jobflow/main)](https://codecov.io/gh/materialsproject/jobflow/) -->
<!-- [![pypi version](https://img.shields.io/pypi/v/jobflow?color=blue)](https://pypi.org/project/jobflow/) -->
<!-- ![supported python versions](https://img.shields.io/pypi/pyversions/jobflow) -->
<!-- [![DOI](https://joss.theoj.org/papers/10.21105/joss.05995/status.svg)](https://doi.org/10.21105/joss.05995) -->

[Documentation](https://battmodels.github.io/asimtools/) | [GitHub](https://github.com/BattModels/asimtools)

# Atomic SIMulation Tools

This package is a lightweight workflow and simulation manager for reproducible
atomistic simulations that can be transferred across environments, calculators
and structures on Unix systems. By using in-built or user-defined asimmodules
and utilities, users can run/build their own simulation recipes and
automagically scale them locally or on slurm based clusters. The core idea is
to separate the dependence of the atomistic potential/calculator and the
simulations steps thereby allowing the same simulation to be run with multiple
calculators/codes and the same calculator to be used for multiple simulation
parameters without altering simulation code. Input and output files follow a
simple consistent file structure and format so that consistent analysis
pipelines can be used across users. For a concrete example of how ASIMTools
achieves this, see the [Developing Custom Asimmodules](https://eeg.engin.umich.edu/asimtools/asimplify.html) page.

## Developer philosophy
The goal of asimtools is to push all the complexity of workflow management,
best-practices, file management etc. into the backend such that the everyday
user only has to handle input files for existing workflows and asimmodule files
for workflows they want to implement.This allows the user to focus on doing the
science and designing experiments. The following are the guiding principles for
how asimtools should work:

- Asimmodules should resemble boilerplate ASE code as much as possible.
- Asimmodules should avoid explicitly depending on a calculator
- Asimmodules should not explicitly depend on the context/environment in which 
  they run
- It should be easy to debug individual asimmodules/parts in large workflows.
  In addition, it should always be possible to debug/resubmit jobs without
  using asimtools.
- Input file structure and format should be standard across all asimmodules. In
  addition all input parameters should match how they would like without
  asimtools i.e. do not provide an API!
- Job progress tracking must be incorporated for easy debugging
- Best practices should be built-in e.g. if multiple jobs of the same slurm
  context are submitted simulataneously, it must be a job array.

## Philosophy on User Experience
The philosophy is to build "simulations" using building blocks of asimmodules.
Asimmodules are nothing but Python functions that return a dictionary, anything
can be done in the function code. These asimmodules can be as
complicated/efficient as you make them using any external packages you want and
can be optimized with time but can still be run within the framework. This
allows a test friendly way to transition from say a tutorial on the
ASE/pymatgen website to an asimtools asimmodule. So while complicated wrappers
are discouraged, they would still work as long as the asimmodule works. The
benefit out of the box is that you can make your asimmodule independent of
calculator or input structures and submit them easily.

We also aim to provide a standard set of robust and efficient simulation
protocols as we develop. You can see all the implemented workflows provided
with the package in the examples directory and modify them to your application.
If you have suggestions for improvements in methodology or code, please bring
up an issue on github!

## Getting Started

These instructions will give you a copy of the project up and running.

### Installing ASIMTools

If you prefer to use a conda a environment, you can create and activate one
using: 
```
conda create -n asimtools python=3.9
conda activate asimtools
```

Then you can install asimtools either using pip or by cloning the source code
from github. Note that the cloned version might have some minor bug fixes that
are not included in the official PyPI release. It is also easier to go through
the examples if you clone the repository.

To install the latest asimtools release from PyPI, you can simply use

```
pip install asimtools
```

To install from source use

```
git clone
https://gitlab.com/ase/ase.git && cd ase pip install . cd ../

git clone https://github.com/BattModels/asimtools.git
cd asimtools
pip install .
```

We recommend you use the latest version of ASE since the ones on PyPI and
conda are quite outdated.

```
git clone https://gitlab.com/ase/ase.git
cd ase
pip install .
```

You can optionally install the dependencies used for development using

```
pip install ".[dev]"
```
Or install some popular Universal Machine Learning Interactomic Potentials
(matgl, mace-torh, chgnet) using

```
pip install ".[mlip]"
```
Making sure the correct versions and dependencies are installed correctly is
probably more stable if you follow their individual installation instructions.

Other individual calculators may need external packages for running them. For
example if you want to use Quantum Espresso or CASTEP, you will have to install
them. Similarly some asimmodules e.g. `lammps.py` might also need external
packages to be used. It is up to the user to make sure those are installed.

You will also need to setup some environment variables, these variables point
to global `env_input.yaml` and `calc_input.yaml` files with your favorite
configurations since these are commonly shared among simulations. You can also
directly specify them when running `asim-execute` (See `asim-execute -h`). You
can also provide a directory for ASIMTools to search for your custom
asimmodules. Examples of these files can be found in the examples.

Add the following to your `.bashrc`
```
export ASIMTOOLS_ENV_INPUT=/path/to/my/global/env_input.yaml
export ASIMTOOLS_CALC_INPUT=/path/to/my/global/calc_input.yaml
export ASIMTOOLS_ASIMMODULE_DIR=/path/to/my/asimmodule/dir
```

## Running the tests

To run tests for the workflow tools, from the tests directory, call:

    pytest

To run the test suite on a component `component.py` , call:

    pytest test_component.py

To run all tests for the provided asimmodules, `cd` into the examples directory
and call:

    source run_all.sh

Or you can run a test for each individual example in its directory using:

    source run.sh

If no errors are reported, the tests have passed. These tests simply test the
functionality of the code, not the scientific validity of the simulations!

To run tests for the provided asimmodules using slurm, you can similarly run
the bash scripts ending with `_slurm.sh`. This will submit a number of jobs,
none longer than 1 minute.
    
## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code
of conduct, and the process for submitting pull requests to us.

## Authors

  - **Keith Phuthi**
    [mkphuthi](https://github.com/mkphuthi)
  - **Emil Annevelink**
    [emilannevelink](https://github.com/emilannevelink)

See also the list of
[contributors](https://github.com/BattModels/asimtools.git/contributors)
who participated in this project.

## License

This project is licensed under the MIT License

<!-- ## Acknowledgments

  - Hat tip to anyone whose code is used -->
