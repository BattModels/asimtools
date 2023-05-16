# Atomic SIMulation Tools

This package is for optimizing and standardizing production runs for atomistic simulations. By using in-built or user-defined scripts,
users can run their own simulation recipes and scale them on slurm based clusters. The core idea is to separate the dependence of the
atomistic potential/calculator and the simulations steps thereby allowing the same simulation to be run with multiple calculators and
the same calculator to be used for multiple simulations. Input and output files must follow a strict format so that consistent 
analysis pipelines can be used across users

## Getting Started

These instructions will give you a copy of the project up and running on
your local machine for development and testing purposes.

### Installing

```
conda install -n asimtools python=3.9
conda activate asimtools

conda install ase -c conda-forge

pip install .
```

This installs the base package. Individual calculators may need external packages
for loading those calculators. It is up to the user to make sure those are installed.

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

    pytest test_*.py
    
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
