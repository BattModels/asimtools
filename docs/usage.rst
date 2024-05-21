Running an Existing Asimmodule
==============================

This section will guide you through running an already existing asimmodule
whether it be one of the asimmodules provided in the main repository, from the
asimmodule repository or a custom asimmodule you have written yourself.

To run a simulation, all you need to do is run either of the following
commands:

.. code-block:: console

    asim-execute sim_input.yaml -c calc_input.yaml -e env_input.yaml

Providing ``calc_input.yaml`` or ``env_input.yaml`` is optional. If not
provided, the globally configured files will be used. See :ref:`envinput`. This
command will automatically run the specified simulation in the correct
directory and environment. 

In rare cases, you can use the following command:

.. code-block:: console

    asim-run sim_input.yaml -c calc_input.yaml

This command will run the simulation in the current directory and environment.
For most cases, you will only ever use ``asim-execute``. The differences
between ``asim-execute`` and ``asim-run`` are explained in
:ref:`asimrunvsasimexecute`.

.. _inputs:

In summary, the steps you need to take to start using any of the in-built
asimmodules, described in detail below, are the following:

1. Install ASIMTools in you python environment or add it to your
   ``PYTHONPATH``.
2. Setup your global ``env_input.yaml`` and ``calc_input.yaml`` and set the
   environment variables pointing to them.
3. Write a ``sim_input.yaml`` based on the examples provided in the repository.
   Do not hesitate to submit an issue if you are confused as we are still in a
   testing phase
4. asim-execute!

Input files
***********

.. _envinput:

env_input.yaml
--------------
One can provide an ``env_input.yaml`` file that details the kind
of environments in which asimmodules can be run. This file can be provided with
the asim-execute command using the ``-e`` flag or set globally. An example of
an env_input.yaml file is given below

.. code-block:: yaml

    # template
    env_id:
      mode:
        use_slurm: true
        interactive: false
        run_prefix: ...
        run_suffix: ...
      slurm:
        flags: [flag1, flag2, ...]
        precommands: [precommand1, precommand2, ...]
        postcommands: [postcommand1, postcommand2, ...]
    
    # Concrete examples below
    inline: # Run the asimmodule directly in the console
      mode:
        use_slurm: false
        interactive: true

    batch_job: # Submit a batch job using slurm with 2 tasks
      mode:
        use_slurm: true
        interactive: false
      slurm: 
        flags:
            - -n 2
            - --mem-per-cpu=2G
        precommands:
            - source ~/.bashrc
            - conda activate asimtools
        postcommands:
            - conda deactivate asimtools
    # Submit an interactive job using slurm, you can use a dictionary
    # for the flags
    interactive_job:
      mode:
        use_slurm: true
        interactive: true
      slurm:
        flags:
          -n: 2
          --gres: gpu:2
        precommands:
          - module load lammps

The highest level key is the ``env_id`` which is the one specified in the
``sim_input.yaml``. An ``env_input.yaml`` can have any number of ``env_id`` s.
That way you can specify one global file if you use the same environments
repeatedly. You can configure a global config file by setting
the environment variable.

.. code-block:: console

    export ASIMTOOLS_ENV_INPUT=/path/to/my/global/env_input.yaml

If you do not provide an ``env_input.yaml`` and there is no file called
``env_input.yaml`` in the work directory, ASIMTools will look for the 
``env_id`` in the global file.

The parameters, required, shown in the template section are described below

- **env_id**: (str) unique key for identifying the environment, ``env_id`` in
  ``sim_input.yaml`` must match one of the ``env_id`` s defined in the 
  ``env_input.yaml`` being used.
- **env_id.mode.use_slurm**: (bool) whether or not to request a slurm
  allocation to run the asimmodule
- **env_id.mode.interactive**: (bool) whether or not to run the asimmodule
  directly in the terminal (using ``salloc``) or to submit a batch job (using
  ``sbatch``).
- **env_id.mode.run_prefix**: (str) string to append before running the
  asimmodule e.g. if ``run_prefix=mpirun`` the asimmodule will be invoked with
  the equivalent of ``mpirun python my_asimmodule.py``. ``run_prefix`` in
  ``env_input.yaml`` is always prepended before the one provided by
  ``calc_input.yaml``.
- **env_id.mode.run_suffix**: (str) string to append after running the
  asimmodule e.g. if ``run_suffix: ' &> out.txt'`` is provided, the asimmodule
  will be invoked with the equivalent of ``python my_asimmodule.py &>
  out.txt``. ``run_suffix`` in ``env_input.yaml`` is always appended after the
  one provided by ``calc_input.yaml``.
- **env_id.slurm.flags**: (list/dict, optional) The slurm flags for the
  allocation as a list of flags e.g. ``[-n 4, -N 1]``. One can also specify a
  dictionary e.g. ``'{-n': 4, '-N': 1, '--mem':2G}``
- **env_id.slurm.precommands**: (list, optional) Commands to be run/added to
  the job script before running the asimmodule. A common use case is loading
  a module or activating an environment.
- **env_id.slurm.postcommands**: (list, optional) Commands to be run/added to
  the job asimmodule after running the asimmodule. e.g. for file cleanup or
  moving files after the job is complete.

.. _calcinput:

calc_input.yaml
---------------
The ``calc_input.yaml`` is used to configure an ASE calculator. As
above, a global configuration file can be set using

.. code-block:: console

    export ASIMTOOLS_CALC_INPUT=/path/to/my/global/calc_input.yaml

or provided to asim-execute at run time. Note that if you launch a chained
workflow with ``asim-run`` instead of ``asim-execute``, asimmodules farther
down the chain will use the global ``calc_input.yaml``, so always use
``asim-execute``

.. code-block:: yaml

  # Template
  calc_id:
    name: ...
    module: ...
    precommands: [precommand1, ...]
    postcommands: [postcommand1, ...]
    run_prefix: ...
    run_suffix: ...
    args:
      arg1: value_1
      ...

  # Concrete examples
  # Here is a simple LJ potential from ASE
  lj: 
    name: LennardJones
    module: ase.calculators.lj
    args:
      sigma: 3.54
      epsilon: 0.00802236

  # GPAW needs a run_prefix to work in parallel using mpirun
  gpaw:
    name: GPAW
    module: gpaw.calculator
    run_prefix: mpirun 
    args:
      kpts: [2,2,2]
      h: 0.1
      xc: PBE
      txt: gpaw_output.txt

  # You can install a universal potential like MACE and define it as well, see
  # asimtools/calculators.py for implemented external calculators. Submit an
  # issue if you want one to be implemented.
  MACE32-medium:
    name: MACE
    args:
      model: medium
      use_device: cuda

The parameters for the calculators provided directly in ASE are specified under
the assumption that the calculator will be initiated as follows:

.. code-block::

    from module import name
    calc = name(**args)

This works for all calculators defined in ASE v3.22 and below, for newer
versions of ASE, you might need to use the versions that use profiles e.g. use
``name: EspressoProfile`` not ``name: Espresso`` until those become stable in
ASE. For externally defined calculators, you can submit an issue and we will
implement it. For example, calculators for NequIP, Deep Potential, MACE, CHGNet
and M3GNet force fields are implemented.

- **calc_id**: (str) unique key for identifying the calculator, ``calc_id`` in
  ``sim_input.yaml`` must match one of the ``calc_id`` s defined in the
  provided ``calc_input.yaml``
- **calc_id.name**: (str) Either the name of the class or the reference to one
  of the provided external calculators. 
- **calc_id.module**: (str) The module from which the calculator class is
  imported. e.g. if ``name=LennardJones`` and ``module=ase.calculators.lj``,
  then the calculator object is imported as ``from ase.calculators.lj import
  LennardJones``. This works if the calculator is available in ASE or follows
  ASE format for initialization such as GPAW. Any other ASE calculator will
  need to have the instantiation defined in :ref:calculators.py
- **calc_id.mode.run_prefix**: (str) string to append before running the
  asimmodule e.g. if ``run_prefix=mpirun`` the asimmodule will be invoked with
  the equivalent of ``mpirun python my_asimmodule.py``. ``run_prefix`` in
  ``env_input.yaml`` is always prepended before the one provided by
  ``calc_input.yaml``.
- **calc_id.mode.run_suffix**: (str) string to append after running the
  asimmodule e.g. if ``run_postfix=' &> out.txt'`` the asimmodule will be
  invoked with the equivalent of ``python my_asimmodule.py &> out.txt``.
  ``run_postfix`` in ``env_input.yaml`` is always appended after the one
  provided by ``calc_input.yaml``.
- **calc_id.precommands**: (list, optional) Commands to be run/added to the job
  asimmodule before running the asimmodule. A common use case is loading a
  module or activating an environment
- **calc_id.postcommands**: (list, optional) Commands to be run/added to the
  job asimmodule after running the asimmodule. e.g. cleaning up bulky tmp or
  wavefunction files
- **calc_id.args**: (dict) key-value pairs to be passed as arguments for the
  initialization of the calculator class. e.g. if the class is LennardJones,
  the arguments are passed as ``calc = LennardJones(**{'sigma':3.2,
  'epsilon':3})``

.. _siminput:

sim_input.yaml
--------------

The minimal requirement to run an asimmodule is to provide a ``sim_input.yaml``
file. An example of a ``sim_input.yaml`` is shown below:

.. code-block:: yaml

    asimmodule: singlepoint 
    env_id: inline
    overwrite: false
    submit: true
    workdir: results
    precommands:
        - export MY_ENV_VAR=3
    args:
        arg1: value_1
        arg2: value_2
        ... 

The parameters are:

- **asimmodule**: (str) name of core asimmodule or /path/to/my/asimmodule.py.
  Core asimmodules defined in the asimmodules directory can be simply referred
  to using Python dot notation. E.g. to specify the
  :func:`asimtools.asimmodules.workflows.sim_array` asimmodule, you would
  specify `workflows.sim_array`. Any other asimmodule should be specified as
  either a full path or a path relative to ``ASIMTOOLS_ASIMMODULE_DIR``
  variable to a python file. E.g. ``my_asimmodules/asim_ple.py``
- **env_id**: (str, optional) Environment/context in which to run asimmodule
  configured in env_input.yaml, defaults to running in the current console
- **overwrite**: (bool, optional) (bool) whether or not to overwrite work
  directories if they exist, defaults to false 
- **submit**: (bool, optional) whether to run the asimmodule. If set to false
  it will just write the input files which is very useful for testing before
  submitting large workflows. You can go in and test one example before
  resubmitting with ``submit=True``, defaults to true 
- **workdir**: (str, optional) The directory in which the asimmodule will be
  run, `asim-execute` will create the directory whereas `asim-run` ignores this
  parameter entirely, defaults to './results'
- **precommands**: (list, optional) a list of commands to run in the console
  before running the asimmodule, defaults to empty list
- **postcommands**: (list, optional) a list of commands to run in the console
  after running the asimmodule, defaults to empty list
- **args**: (dict) The arguments of the function being called in the asimmodule
  as key-value pairs. These are specific to the asimmodule being run.

All ASIMTools generated files are named ``sim_input.yaml`` but you can name
user defined files as whatever you like

.. _specifyingimages:

Specifying Images/Atoms
-----------------------

One of the most useful applications of ASIMTools is the unification of methods
for setting up ASE atoms objects using the same interface. If an asimmodule
requires a single or multiple atoms objects as input, they are provided as
either an ``image`` dictionary for a single Atoms object or ``images`` for a
list of Atoms objects as part of the ``args`` section. Below are the different
ways to get an atoms object. You can also download images from The Materials
Project and for some cases generate them using Pymatgen.

For a detailed description of the API and examples, see
:func:`asimtools.utils.get_atoms`

.. code-block:: yaml

  # Reading a specific image from a structure file using ase.io.read
  image:
    image_file: /path/to/my/ASE-readable/image/file.xyz
    # Optional keyword argument passed to ase.io.read
    index: 3
  
  # Building a bulk crystal using ase.build.bulk
  image:
    builder: bulk
    # Optional keyword arguments passed to the builder, must match ASE exactly
    name: Li
    crystalstructure: bcc
    a: 4.3
    cubic: True

  # Building a surface using ase.build.fcc100
  image:
    builder: fcc100
    # Optional keyword arguments passed to the builder, must match ASE exactly
    symbol: Fe
    vacuum: 8
    periodic: False

  # Building a 3x3x3 supercell of Ar using ase.build.bulk then
  # Atoms.repeat(repeat) and then applying Atoms.rattle(stdev=rattle_stdev)
  image:
    name: Ar
    repeat: [3,3,3]
    rattle_stdev: 0.01

  # You can even supply an atoms object directly so that the interface is
  # universal. This is most useful in the asimmodule code itself.
  image:
    atoms: Atoms

  # An example downloading a structure from Materials Project using your own
  # USER_API_KEY
  image:
    mp_id: 'mp-14'
    interface: pymatgen
    user_api_key: "USER_API_KEY"
    conventional_unit_cell': true


Similarly, if the asimmodule requires multiple image inputs, there exists a
universal interface. The keyword is usually specified as ``images``. This is
especially useful for distributing simulations across multiple structures or
reading structures from multiple previous simulations, even in different
directories.

For a detailed description of the API, see :func:`asimtools.utils.get_images`

.. code-block:: yaml

  # Reading specific images from a structure file using ase.io.read
  images:
    image_file: /path/to/my/ASE-readable/image/file.xyz
    # Optional keyword arguments passed to ase.io.read
    index: '3:8'
    format: extxyz
  
  # You can read all files matching a certain pattern using a wildcard
  images:
    pattern: /path/to/my/structure/files/*.cif
    # Optional keyword argument passed to ase.io.read
    index: -1

  # You can read all files matching certain patterns using a wildcard
  images:
    patterns: 
    - /path/to/my/structure/files/*.cif
    - /path/to/my/other/structure/files/*.cfg
  
  # You can even supply a list of atoms objects directly so that the interface
  # is universal. This is most useful in the asimmodule code itself.
  images:
    images: [Atoms1, Atoms2, ...]

.. _asimrunvsasimexecute:

Usage of asim-execute and asim-run  
**********************************
The major difference between ``asim-execute`` and ``asim-run`` is that,
``asim-execute`` takes into account the ``workdir`` and the ``env_id``.
``asim-run`` will run the asimmodule in the current directory and in the
current console. In fact, ``asim-execute`` will create the ``workdir`` and then
run ``asim-run`` in the correct environment/batch job. You can always for
example, request a slurm allocation, go to the directory where you want the
asimmodule to be run and call ``asim-run`` from there if you would like more
control or to debug. If you want verbose logs for debugging, you can run with
the ``-d`` or ``--debug`` flag.

.. _outputs:

Output files
************
A job or asimmodule run through ASIMTools will always produce a standard set of
output files in addition to whatever outputs the asimmodule produces. In
particular the most important outputs are the ``output.yaml`` and the
``job.log`` file. 

#. \``output.yaml`` contains the status of the job being run in the current
   directory which can be one of ``clean, started, complete, failed, discard``.
   The statuses are self-explanatory, the ``discard`` status is never written
   by ASIMTools but a user can edit an ``output.yaml`` file and change it's
   status to ``discard`` to tell ASIMTools to ignore that job in any workflows.
   This is common for example if you launch multiple jobs and one of them fails
   irredemably. Deleting the directory for that job is also ok if nothing
   depends on it downstream. Importantly, any results returned by the function
   defined in the asimmodule are found in ``output.yaml``. Asimmodule functions
   should always return a dictionary of only primitive types for this purpose.

   An example of an ``output.yaml`` file is shown below.

.. code-block:: yaml

  # Successful output for singlepoint asimmodule
  end_time: 2023-08-28 21:50:51.368025
  energy: 13.77302319846367  #This was added by the scinglepoint asimmodule
  files:
    image: image_output.xyz
  job_ids: '372919'
  start_time: 2023-08-28 21:50:46.188300
  status: complete

  # Failed output
  start_time: 14:29:55, 10/06/23
  status: failed

#. ``job.log`` captures the logged output of ``asim-run`` or asimmodules
   that use logging. It is extremely useful for debugging as following the logs
   starting from the base directory will usually lead you to the correct
   traceback that caused the failure.

#. ``stderr.txt`` captures errors and backtraces from running asimmodules. This
   is usually the most informative file for debugging. You can be directed to
   the correct one by noting errors in ``job.log`` files.

#. ``stdout.txt`` captures any stdout from running asimmodules. It is mostly a
   safety measure for catching anything that prints to stdout and rarely has
   useful information unless you write an asimmodule that uses ``print``
   statements. In batch jobs, this output this goes to the slurm job output.

#. ``input_image.xyz`` and ``input_images.xyz`` capture the images input into
   the asimmodule. This makes sure there is a concrete artifact for the
   structure used by the asimmodule for the purposes of visualization and
   debugging. They are always in ``extxyz`` format as a flexible standard
   format

#. ``slurm*`` are slurm job files which can be named according to flags
   specified in ``env_input.yaml`` otherwise are named
   ``slurm_stdout.id-%a_j%A`` or ``slurm_stderr.id-%a_j%A`` after job and array
   IDs

.. _restarting:

Checking job status and Restarting failed jobs
**********************************************
To check the status of jobs, even complicated chains and distributed jobs, we
provide the ``asim-check`` utility which can be run using:

.. code-block:: console

  asim-check /path/to/sim_input.yaml

This will print the job tree, including statuses and work directories of the
jobs whose root directory is specified as ``workdir`` in ``sim_input.yaml``.

In many cases, there may be mistakes in one of your configuration files leading
to a failed workflow. In these cases there are a few ways you could resolve the
issue:

* Delete the work directory and restart the workflow. This is why it is
  recommended that the base ``sim_input.yaml`` has ``workdir`` set to a new
  directory that only has the results of the workflow.
* Modify the ASIMTools generated ``sim_input.yaml`` to fix the problem. If
  there are downstream ``sim_input.yaml`` files in a chain, they will have to
  be deleted or set ``overwrite=True``. Deleting is recommended for safety
  purposes.
* Simply rerun ``asim-execute``. This will rerun the jobs, skipping any jobs
  with a status of ``complete`` or ``discard``. Note that error files are not
  deleted so you will have to clear those manually. Use this with caution!

Importing functions from asimmodules and the API
************************************************

Because asimmodules contain what are merely Python functions, you can always
import them and use them in any other code for example, you can import
:func:`asimtools.asimmodules.singlepoint` and use it as below.

.. code-block:: python
  
  from asimtools.asimmodules.singlepoint import singlepoint

  results = singlepoint(image={'name': 'Ar'}, calc_id='lj')
  print(results)

You can also use the utils and tools e.g. to load a calculator using just a
``calc_id``

.. code-block:: python

  from asimtools.calculators import load_calc

  calc = load_calc('lj')
