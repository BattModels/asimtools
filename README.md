# Atomic SIMulation Tools

This package is for optimizing and standardizing production runs for atomistic simulations. By using in-built or user-defined scripts,
users can run their own simulation recipes and scale them on slurm based clusters. The core idea is to separate the dependence of the
atomistic potential/calculator and the simulations steps thereby allowing the same simulation to be run with multiple calculators and
the same calculator to be used for multiple simulations. Input and output files must follow a strict format so that consistent 
analysis pipelines can be used across users

## Getting Started

These instructions will give you a copy of the project up and running on
your local machine for development and testing purposes.

### Prerequisites

Requirements for the software and other tools to build, test and push 
- [ASE](https://wiki.fysik.dtu.dk/ase/index.html)
- [PyMatgen](https://pymatgen.org/)

In addition to the universal requirements, individual calculators may need external packages
for loading those calculators. It is up to the user to make sure those are installed.

### Installing

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

    pytest component.py

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