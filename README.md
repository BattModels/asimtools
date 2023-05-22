# Atomic SIMulation Tools

This package is for optimizing and standardizing production runs for atomistic simulations. By using in-built or user-defined scripts,
users can run their own simulation recipes and scale them on slurm based clusters. The core idea is to separate the dependence of the
atomistic potential/calculator and the simulations steps thereby allowing the same simulation to be run with multiple calculators and
the same calculator to be used for multiple simulations. Input and output files must follow a strict format so that consistent 
analysis pipelines can be used across users

## Developer philosophy
The main goal of asimtools is to push all the complexity of workflow management, best-practices, file management etc. into the backend such that the everyday user only has to handle input files for existing workflows and script files for workflows they want to implement. The following are the guiding principles for how asimtools should work:

- Scripts should resemble boilerplate ASE code as much as possible.
- Scripts should not explicitly depend on a calculator
- Scripts should not explicitly depend on the context in which they run
- It should be easy to debug individual scripts/parts in large workflows. In addition, it should always be possible to debug/resubmit jobs without using asimtools.
- input file structure and format should be staandard across all scripts. In addition all input parameters should match how they would like without asimtools i.e. do not provide an API! This means extensive use of importlib.
- output file structure should be predictable across all scripts
- Job progress tracking must be incorporated
- Best practices should be built-in e.g. if multiple jobs of the same context are submitted simulataneously, it must be a job array

## Getting Started

These instructions will give you a copy of the project up and running on
your local machine for development and testing purposes.

### Installing
You can install asimtools in a new conda environment using:
```
conda create -n asimtools python=3.9
conda activate asimtools

conda install ase -c conda-forge

git clone https://github.com/BattModels/asimtools.git
cd asimtools
pip install -e .
```

This installs the base package in developer mode so that you do not have to
reinstall every time you make changes.

Individual calculators may need external packages for loading those calculators. It is up to the user to make sure those are installed.

## Examples
TODO: A step by step series of examples that tell you how to get a development
environment running

Say what the step will be

    Give the example

And repeat

    until finished

End with an example of getting some data out of the system or using it
for a little demo


## Running the tests

To run all tests from the tests directory, call:

    pytest

To run the test suite on a component `component.py` , call:

    pytest test_component.py
    
## Basic example

Simulations are run by providing a `*calc_input.yaml` and `*sim_input.yaml` file which specify 
the calculator (and the environment it runs in) and the simulation parameters which are specific 
to the simulation being run. The recommended method for calling scripts is to use

```
asim-run *calc_input.yaml *sim_input.yaml
```

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code
of conduct, and the process for submitting pull requests to us.

## Authors

  - **Keith Phuthi**
    [mkphuthi](https://github.com/mkphuthi)

See also the list of
[contributors](https://github.com/BattModels/asimtools.git/contributors)
who participated in this project.

## License

This project is licensed under the ABC License

## Acknowledgments

  - Hat tip to anyone whose code is used
