Running an existing asimmodule
==========================

This section will guide you through running an already existing asimmodule whether
it be one of the asimmodules provided in the main repository, from the asimmodule
repository or a custom asimmodule you have written yourself.

To run a simulation, all you need to do is run either of the following
commands:

.. code-block:: console

    asim-execute sim_input.yaml -c calc_input.yaml -e env_input.yaml

or

.. code-block:: console

    asim-run sim_input.yaml -c calc_input.yaml

In either case, providing ``calc_input.yaml`` or ``env_input.yaml`` is optional. If
not provided, the globally configured files will be used. See :ref:`envinput`. For most
cases, you will only ever use ``asim-execute``. The differences between
``asim-execute`` and ``asim-run`` are explained in :ref:`asimrunvsasimexecute`.

.. _inputs:

Input files
***********

.. _siminput:

sim_input.yaml
--------------

The minimal requirement to run a asimmodule is to provide a ``sim_input.yaml`` file.
An example of a ``sim_input.yaml`` is shown below:

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

- **asimmodule**: (str) name of core asimmodule or /path/to/my/asimmodule.py 
- **env_id**: (str, optional) Environment/context in which to run asimmodule
  configured in env_input.yaml, defaults to running in the current console
- **overwrite**: (bool, optional) (bool) whether or not to overwrite work
  directories if they exist, defaults to false 
- **submit**: (bool, optional) whether or not to run the asimmodule or instead just
  write the files to run the asimmodule without running the asimmodule, defaults to true 
- **workdir**: (str, optional) The directory in which the asimmodule will be run,
  `asim-execute` will create the directory whereas `asim-run` ignores this
  parameter, defaults to './'
- **precommands**: (list, optional) a list of commands to run in the console
  before running the asimmodule, defaults to empty list
- **postcommands**: (list, optional) a list of commands to run in the console
  after running the asimmodule, defaults to empty list
- **args**: (dict) The arguments of the function being called in the asimmodule as
  key-value pairs. These are specific to the asimmodule being run.

All ASIMTools generated files are named ``sim_input.yaml`` but you can name
user defined files as whatever you like

.. _envinput:

env_input.yaml
--------------
Additionally, one can provide an ``env_input.yaml`` file that details the kind of
environments in which asimmodules can be run. This file can be provided with the
asim-execute command using the ``-e`` flag or set globally. An example of an
env_input.yaml file is given below

.. code-block:: yaml

    # template
    env_id:
      mode:
        use_slurm: true
        interactive: false
        run_prefix: ...
        run_postfix: ...
      slurm:
        flags: [flag1, flag2, ...]
        precommands: [precommand1, precommand2, ...]
        postcommands: [postcommand1, postcommand2, ...]
    
    # Concrete examples below
    inline: # Run the asimmodule directly in the commandline
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
        precommands:
            - source ~/.bashrc
            - conda activate asimtools
        postcommands:
            - conda deactivate asimtools

    interactive_job: # Submit an interactive job using slurm
      mode:
        use_slurm: true
        interactive: true
      slurm:
        flags:
          - -n 2
          - --gres:gpu=2
        precommands:
          - module load lammps

The highest level key is the ``env_id`` which is specified in the
``sim_input.yaml``. An ``env_input.yaml`` can have any number of ``env_id`` s. That
way you can specify one global file if you use the same environments
repeatedly. In particular, you can configure a global config file by setting
the environment variable.

.. code-block:: console

    export ASIMTOOLS_ENV_INPUT=/path/to/my/global/env_input.yaml

The parameters, required, shown in the template section are  are described below

- **env_id**: (str) unique key for identifying the environment, ``env_id`` in
  ``sim_input.yaml`` must match one of the ``env_id`` s defined in the provided
  ``env_input.yaml``
- **env_id.mode.use_slurm**: (bool) whether or not to request a slurm
  allocation to run the asimmodule
- **env_id.mode.interactive**: (bool) whether or not to request a slurm
  allocation to run the asimmodule directly in the terminal (using ``salloc``) or
  to submita batch job (using ``sbatch``)
- **env_id.mode.run_prefix**: (str) string to append before running the asimmodule
  e.g. if ``run_prefix=mpirun`` the asimmodule will be inkoked with the equivalent
  of ``mpirun python my_asimmodule.py``. ``run_prefix`` in ``env_input.yaml`` is
  always prepended before the one provided by ``calc_input.yaml``.
- **env_id.mode.run_postfix**: (str) string to append after running the asimmodule
  e.g. if ``run_postfix=' &> out.txt'`` the asimmodule will be inkoked with the
  equivalent of ``python my_asimmodule.py &> out.txt``. ``run_postfix`` in
  ``env_input.yaml`` is always appended after the one provided by
  ``calc_input.yaml``.
- **env_id.slurm.flags**: (list, optional) The slurm flags for the allocation
  as a list of flags e.g. ``[-n 4, -N 1]``
- **env_id.slurm.precommands**: (list, optional) Commands to be run/added to
  the job asimmodule before running the asimmodule. A common use cas is loading a
  module or activating an environment
- **env_id.slurm.postcommands**: (list, optional) Commands to be run/added to
  the job asimmodule after running the asimmodule.

.. _calcinput:

calc_input.yaml
---------------
Lastly the ``calc_input.yaml`` is used to configure an ASE calculator. As
above, a global configuration file can be set using

.. code-block:: console

    export ASIMTOOLS_CALC_INPUT=/path/to/my/global/calc_input.yaml

or provided to asim-execute at run time. Note that if you launch a chained
workflow with ``asim-run`` instead of ``asim-execute``, asimmodules farther down
the chain will use the global ``calc_input.yaml``


.. code-block:: yaml

  # Template
  calc_id:
    name: ...
    module: ...
    precommands: [precommand1, ...]
    postcommands: [postcommand1, ...]
    run_prefix: ...
    run_postfix: ...
    args:
      arg1: value_1
      ...

  # Concrete examples
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


- **calc_id**: (str) unique key for identifying the calculator, ``calc_id`` in
  ``sim_input.yaml`` must match one of the ``calc_id`` s defined in the
  provided ``calc_input.yaml``
- **calc_id.name**: (str) Either the name of the class or the reference to one
  of the provided BBBB external calculators. 
- **calc_id.module**: (str) The module from which the calculator class is
  imported. e.g. if ``name=LennardJones`` and ``module=ase.calculators.lj``,
  then the calculator object is imported as ``from ase.calculators.lj import
  LennardJones``. This works if the calculator is available in ASE or follows
  ASE format for initialization such as GPAW. Any other ASE calculator will
  need to have the instantiation defined in :ref:calculators.py
- **calc_id.mode.run_prefix**: (str) string to append before running the asimmodule
  e.g. if ``run_prefix=mpirun`` the asimmodule will be inkoked with the equivalent
  of ``mpirun python my_asimmodule.py``. ``run_prefix`` in ``env_input.yaml`` is
  always prepended before the one provided by ``calc_input.yaml``.
- **calc_id.mode.run_postfix**: (str) string to append after running the asimmodule
  e.g. if ``run_postfix=' &> out.txt'`` the asimmodule will be inkoked with the
  equivalent of ``python my_asimmodule.py &> out.txt``. ``run_postfix`` in
  ``env_input.yaml`` is always appended after the one provided by
  ``calc_input.yaml``.
- **calc_id.precommands**: (list, optional) Commands to be run/added to the job
  asimmodule before running the asimmodule. A common use cas is loading a module or
  activating an environment
- **calc_id.postcommands**: (list, optional) Commands to be run/added to the
  job asimmodule after running the asimmodule.
- **calc_id.args**: (dict) key-value pairs to be passed as arguments for the
  initialization of the calculator class. e.g. if the class is LennardJones,
  the arguments are passed as ``calc = LennardJones(**{'sigma':3.2,
  'epsilon':3})``

.. _specifyingimages:

Specifying Images/Atoms
-----------------------

One of the most useful applications of ASIMTools is the unification of methods
for setting up ASE atoms objects using the same interface. If a asimmodule requires
a single or multiple atoms objects as input, they are provided as either an
``image`` dictionary for a single Atoms object or ``images`` for a list of
Atoms objects as part of the ``args`` section. Below are all the different ways
to get an atoms object. Downloading images from MaterialsProject and Generating
them from Pymatgen will implemented in future.

For a detailed deasimmoduleion of the API, see :func:`asimtools.utils.get_atoms`

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

Similarly, if the asimmodule requires multiple image inputs, there exists a
universal interface. The keyword is strictly specified as ``images`` This is especially useful for distributing simulations
across multiple structures or reading structures from multiple previous
simulations.

For a detailed deasimmoduleion of the API, see :func:`asimtools.utils.get_images`

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

  # You can read all files matching a certain pattern using a wildcard
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
``asim-run`` will the asimmodule in the current directory and in the current
console. In fact, ``asim-execute`` will create the ``workdir`` and then run
``asim-run`` in the correct environment/batch job. You can always for example,
request a slurm allocation, go to the directory where you want the asimmodule to be
run and call ``asim-run`` from there if you would like more control or to debug.
If you want verbose logs for debugging, you can run with the ``-d`` or 
``--debug`` flag.

.. _outputs:

Output files
****************
A job or asimmodule run through ASIMTools will always produce a standard set of
output files in addition to whatever outputs the asimmodule produces. In particular
the most  important output are the ``output.yaml`` and the ``job.log`` file. 

#. \``output.yaml`` contains the status of the job being run in the current
   directory which can be one of ``clean, started, complete, failed, discard``.
   The statuses are self-explanatory, the ``discard`` status is never written
   by ASIMTools but a user can edit an ``output.yaml`` file and change it's
   status to ``discard`` to tell ASIMTools to ignore that job in any workflows.
   This is common for example if you launch multiple jobs and one of them fails
   irredemably. Deleting the directory for that job is also ok if nothing
   depends on it downstream. Importantly, any results returned by the function
   defined in the asimmodule are found in ``output.yaml``. asimmodule main functions
   should always return a dictionary for this purpose.

   An example of an ``output.yaml`` file is shown below

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

#. ``stderr.txt`` captures errors and backtraces from running asimmodules. This is
   usually the most informative file for debugging. You can be directed to the
   correct one by noting errors in ``job.log``

#. ``stdout.txt`` captures any stdout from running asimmodules. It is mostly a
   safety measure for catching anything that prints to stdout and rarely has
   useful information unless you write a asimmodule that uses ``print`` statements.

#. ``input_image.xyz`` and ``input_images.xyz`` capture the images input into
   the asimmodule. This makes sure there is a concrete artifact for the structure
   used in the asimmodule for the purposes of visualization and debugging. They are
   always in ``extxyz`` format as a flexible standar format

#. ``slurm*`` are slurm job files which can be named according to flags
   specified in ``env_input.yaml`` otherwise are named ``slurm_[job_id].out``


.. _restarting:

Checking job status and Restarting failed jobs
**************************************************
To check the status of jobs, even complicated chains and distributed jobs, we
provide the ``asim-check`` utility which can be run using:

.. code-block:: console

  asim-check /path/to/sim_input.yaml

This will print the job tree, including statuses and work directories of the
jobs whose root directory is specified as ``workdir`` in ``sim_input.yaml``.

In many cases, there may be mistakes in one of your configuration files leading to
a failed workflow. In these cases there are a few ways you could resolve the
issue:

* Delete the workdirectory and restart the workflow. This is why it is
  recommended that the base ``sim_input.yaml`` has ``workdir`` set to a new
  directory that only has the results of the workflow.
* Modify the ASIMTools generated ``sim_input.yaml`` to fix the problem. If
  there are downstream ``sim_input.yaml`` files in a chain, they will have to
  be deleted or set ``overwrite=True``. Deleting is recommended for safety
  purposes.
* Simply rerun ``asim-execute``. This will rerun the jobs, skipping any jobs
  with a status of ``complete`` or ``discard``. Note that error files are not
  deleted so you will have to clear those manually

Importing functions from asimmodules
********************************

Because asimmodules contain what are merely Python functions, you can always import
them and use them in any other code for example, you can import 
:func:`asimtools.asimmodules.singlepoint` and use it as below.

.. code-block:: python
  
  from asimtools.asimmodules.singlepoint import singlepoint

  results = singlepoint(image={'name': 'Ar'}, calc_id='lj')
  print(results)
