# Atomic SIMulation Tools

This package is a lightweight workflow and simulation manager for reproducible
atomistic simulations that can be transferred across environments, calculators
and structures. By using in-built or user-defined asimmodules and utilities, users
can run/build their own simulation recipes and automagically scale them on
slurm based clusters. The core idea is to separate the dependence of the
atomistic potential/calculator and the simulations steps thereby allowing the
same simulation to be run with multiple calculators and the same calculator to
be used for multiple simulations. Input and output files must follow a strict
format so that consistent analysis pipelines can be used across users

** ASIMTools is going through some breaking changes for version 0.1.0, please use with caution **

## Documentation
Documentation is under construction [here](https://battmodels.github.io/asimtools/)

## Developer philosophy
The goal of asimtools is to push all the complexity of workflow
management, best-practices, file management etc. into the backend such that the
everyday user only has to handle input files for existing workflows and asimmodule
files for workflows they want to implement. The following are the guiding
principles for how asimtools should work:

- asimmodules should resemble boilerplate ASE code as much as possible.
- asimmodules should not explicitly depend on a calculator
- asimmodules should not explicitly depend on the context in which they run
- It should be easy to debug individual asimmodules/parts in large workflows. In
  addition, it should always be possible to debug/resubmit jobs without using
  asimtools.
- input file structure and format should be standard across all asimmodules. In
  addition all input parameters should match how they would like without
  asimtools i.e. do not provide an API! This means extensive use of importlib.
- output file structure should be predictable across all asimmodules
- Job progress tracking must be incorporated
- Best practices should be built-in e.g. if multiple jobs of the same slurm
  context are submitted simulataneously, it must be a job array

## Philosophy on User Experience
The philosophy is to build "simulations" using building blocks of asimmodules.
These asimmodules can be as complicated/efficient as you make them using any
external packages you want and can be optimized with time but can still be run
within the framework. This allows a test friendly way to transition from say a
tutorial on the ASE/pymatgen website to an asimtools asimmodule. So while
complicated wrappers are discouraged, they would still work as long as the
asimmodule works. The benefit out of the box is that you can make your asimmodule
independent of calculator or input structures and submit them easily.

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
reinstall if you make changes to the core package.

Individual calculators may need external packages for loading those
calculators. Similarly some asimmodules e.g. `lammps.py` might also need external packages to be used. It is up to the user to make sure those are installed.

You will also need to setup some environment variables, these variables point
to a `env_input.yaml` and `calc_input.yaml` with your favorite configs since
these are commonly shared among simulations. You can also directly specify them
when running `asim-execute` (See `asim-execute -h`). 
Examples of these files can be found in the examples.

Add the following to your `.bashrc`
```
export ASIMTOOLS_ENV_INPUT=/path/to/my/global/env_input.yaml
export ASIMTOOLS_CALC_INPUT=/path/to/my/global/calc_input.yaml
export ASIMTOOLS_asimmodule_DIR=/path/to/my/asimmodule/dir
```
<!-- ## Examples -->
<!-- Below are a few examples on how to run already implemented asimmodules. -->
<!-- 
The first thing to understand is the difference between `asim-execute sim_input.yaml` and 
`asim-run sim_input.yaml`. The latter runs the chosen asimmodule in whatever location 
and environment it happens to be launched from i.e. equivalent to 
`python my_asimmodule.py sim_input.yaml` whereas the latter runs submits the job 
that runs the asimmodule e.g. it will go to the work directory and launch a slurm job 
from there containing `asim-run sim_input.yaml`.  `asim-execute` is essentially 
"environment-aware" and runs the asimmodule in the environment specified by `env_id`
whereas `asim-run` does not use `env_id` at all. -->

<!-- A template is provided for writing a new asimmodule called `template.py`

We provide a standard set of building block asimmodules for simulations/workflows that tend to form components of large simulations. Let us consider some examples. 
To see details of their arguments (`args`), see their docstrings (which don't exist yet :( )

1. The simplest are "unit" asimmodules which do not internally call other asimmodules 
or depend on results from other asimmodules. 
The simplest example is `singlepoint.py` which runs an energy/force/stress calculation
on a single structure. To launch in your environment of choice specified by
`env_id` in `sim_input.yaml` use `asim-execute sim_input.yaml`. See `asim-execute -h`. 
For all asimmodules, you have the ability to specify the `env_input.yaml` and `calc_input.yaml` or 
skip them to use the global files specified in `ASE_CALC_INPUT` and 
`ASE_ENV_INPUT`. Note that many asimmodules don't need `calc_input.yaml` if they 
just do some preprocessing or analysis. Another example is `atom_relax.py`.
*Exercise: Add `cell_relax.py` based on `atom_relax.py`*

2. Another type are "parent" asimmodules which launch multiple "child" asimmodules in parallel.
An example is `image_array.py` which runs the same simulation on multiple
images e.g. you want the energies of all the structures in a database. Note
that you only need to point to the images, the env and the calc in the sim_input,
The asimmodule, which uses a DistributedJob object automatically handles whether to
submit jobs using a slurm array, individual slurm jobs or one after the other
depending on the specified `env_id`. Another example is 
`strong_scaling.strong_scaling.py`. 
*Exercise: Implement `env_array.py` and `calc_array.py` based on `image_array.py` and `strong_scaling.py`.*

3. The last major type of "parent" asimmodule are chained asimmodules. These asimmodules run one 
**UnitJob** after the other based on a specified workflow. The job automatically runs
one job after another or submits slurm jobs with appropriate dependencies. Note
that if one of the unitjobs internally calls another asimmodule such as 
`image_array.py`, then the slurm dependencies will fail but everything else should work.
In this casea, the current workaround is to set `submit=false` for the stages of 
the chain that would fail due to a previous array/chain job. Then rerun `asim-execute`
again with `submit=true`, It should automatically skip the completed steps. See `eos.eos.py` 
and the corresponding example. An example to help understand the genera; chained sim_input is 
`chained`

4. Hybrid jobs combine multiple of these elements e.g. `eos.eos.py` runs a
`preprocessing` asimmodule to prepare the scaled structures followed by a UnitJob that submits
an `image_array` for each scaled structure and finally a unitjob to `postprocess`
the results from the array. It is possible to 
- build hybrid asimmodules in python without using existing asimtools asimmodules
- running a `chained.py` on a `sim_input.yaml` with the correct format that runs asimmodules that are already defined
- to write a asimmodule such that it directly manipulates Job objects. This is the most flexible and robust but ofcourse most complicated
*Exercise: Run the eos calculation without writing any new asimmodules, simply using the asimmodules in the core (image_array) and eos subfolder (preprocessing and postprocessing).* 
The key is being able to to point to the correct files for each step before the calculation is run. -->

## Running the tests

To run all tests from the tests directory, call:

    pytest

To run the test suite on a component `component.py` , call:

    pytest test_component.py
    
<!-- ## Basic example

Simulations are run by providing a `*calc_input.yaml` and `*sim_input.yaml` file which specify 
the calculator (and the environment it runs in) and the simulation parameters which are specific 
to the simulation being run. The recommended method for calling asimmodules is to use

```
asim-run *calc_input.yaml *sim_input.yaml
``` -->

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
