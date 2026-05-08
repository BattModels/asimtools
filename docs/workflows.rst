Using built-in workflow tools
=============================

ASIMTools provides six built-in workflow asimmodules for the most common
simulation patterns. They are all wrappers around the
:class:`~asimtools.job.DistributedJob` and :class:`~asimtools.job.ChainedJob`
classes and can be nested arbitrarily.

.. contents:: Workflows
   :local:
   :depth: 1

----

sim\_array
----------

**Module:** :func:`asimtools.asimmodules.workflows.sim_array.sim_array`

Run the **same asimmodule in parallel** over a sweep of values for a single
argument (or a paired set of arguments). This is the most generally useful
workflow and covers parameter sweeps, convergence tests, and ensemble runs.

When to use
~~~~~~~~~~~

- Sweep a lattice constant, cutoff energy, temperature, etc.
- Vary an integer index (e.g. run the same relaxation on structures 0-99).
- Sweep multiple parameters simultaneously (``secondary_key_sequences``).

Key parameters
~~~~~~~~~~~~~~

+-----------------------------+--------------------------------------------------+
| Parameter                   | Description                                      |
+=============================+==================================================+
| ``template_sim_input``      | Base sim_input copied for every job.             |
+-----------------------------+--------------------------------------------------+
| ``key_sequence``            | Dot-path into the nested dict being swept,       |
|                             | e.g. ``[args, a]`` to change ``args.a``.         |
+-----------------------------+--------------------------------------------------+
| ``array_values``            | Explicit list of values.                         |
+-----------------------------+--------------------------------------------------+
| ``linspace_args``           | ``[start, stop, n]`` — passed to                 |
|                             | :func:`numpy.linspace`.                          |
+-----------------------------+--------------------------------------------------+
| ``arange_args``             | ``[start, stop, step]`` — passed to              |
|                             | :func:`numpy.arange`.                            |
+-----------------------------+--------------------------------------------------+
| ``file_pattern``            | Glob pattern; files matching it become values.   |
+-----------------------------+--------------------------------------------------+
| ``labels``                  | ``"values"`` (default) uses the value itself;    |
|                             | ``"files"`` extracts a label from the file path; |
|                             | a list provides explicit directory names.        |
+-----------------------------+--------------------------------------------------+
| ``placeholder``             | Replace this sentinel string inside the value    |
|                             | rather than the whole value.                     |
+-----------------------------+--------------------------------------------------+
| ``secondary_key_sequences`` | Additional keys to sweep in lock-step.           |
+-----------------------------+--------------------------------------------------+
| ``group_size``              | Pack N jobs into each Slurm array task.          |
+-----------------------------+--------------------------------------------------+

Example — sweep lattice constants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   # sim_input.yaml
   asimmodule: workflows.sim_array
   env_id: batch
   workdir: results/lc_sweep
   args:
     key_sequence: [args, image, a]
     linspace_args: [3.4, 4.2, 9]   # 9 values from 3.4 to 4.2 Å
     template_sim_input:
       asimmodule: singlepoint
       env_id: batch
       args:
         calculator:
           calc_id: my_mlip
         image:
           name: Cu
           crystalstructure: fcc
           a: PLACEHOLDER
         properties: [energy, forces, stress]

Example — sweep over a set of structure files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   asimmodule: workflows.sim_array
   env_id: batch
   workdir: results/structure_sweep
   args:
     key_sequence: [args, image, image_file]
     file_pattern: data/structures/*.xyz
     labels: files             # label taken from the filename
     template_sim_input:
       asimmodule: geometry_optimization.cell_relax
       env_id: batch
       args:
         calculator:
           calc_id: my_mlip
         image:
           image_file: PLACEHOLDER

Example — sweep two parameters simultaneously
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   asimmodule: workflows.sim_array
   env_id: batch
   workdir: results/dual_sweep
   args:
     key_sequence: [args, temperature]
     array_values: [300, 600, 900, 1200]
     secondary_key_sequences:
       - [args, pressure]
     secondary_array_values:
       - [1.0, 2.0, 3.0, 4.0]
     template_sim_input:
       asimmodule: ase_md.ase_md
       env_id: batch
       args:
         calculator:
           calc_id: my_mlip
         image: {name: Fe}
         temperature: 300
         pressure: 1.0

----

image\_array
------------

**Module:** :func:`asimtools.asimmodules.workflows.image_array.image_array`

Run the **same asimmodule in parallel on multiple structures**. The image
input specification for each job is taken from
:func:`asimtools.utils.get_images`, so any combination of file patterns,
explicit lists, or ASE database queries is supported.

When to use
~~~~~~~~~~~

- Relaxation or single-point calculation across a database of structures.
- Benchmark a force field against a set of reference structures.
- Compute phonons for every structure in a dataset.

Key parameters
~~~~~~~~~~~~~~

+-----------------------------+---------------------------------------------------+
| Parameter                   | Description                                       |
+=============================+===================================================+
| ``images``                  | Dict for :func:`~asimtools.utils.get_images`;     |
|                             | can use ``image_file``, ``pattern``,              |
|                             | ``patterns``, or ``images``.                      |
+-----------------------------+---------------------------------------------------+
| ``template_sim_input``      | Base sim_input; image key is injected per job.    |
+-----------------------------+---------------------------------------------------+
| ``key_sequence``            | Where in ``template_sim_input`` to inject the     |
|                             | image path (default ``[args, image,               |
|                             | image_file]``).                                   |
+-----------------------------+---------------------------------------------------+
| ``labels``                  | Source of directory name for each job.            |
+-----------------------------+---------------------------------------------------+
| ``env_ids``                 | Per-image environments (or a single string).      |
+-----------------------------+---------------------------------------------------+

Example — single-point energy on a dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   asimmodule: workflows.image_array
   env_id: batch
   workdir: results/dataset_sp
   args:
     images:
       pattern: data/structures/*.xyz
     template_sim_input:
       asimmodule: singlepoint
       env_id: batch
       args:
         calculator:
           calc_id: my_mlip
         properties: [energy, forces, stress]
     labels: files

Example — relax all structures in an xyz file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   asimmodule: workflows.image_array
   env_id: batch
   workdir: results/relax_all
   args:
     images:
       image_file: data/structures.xyz
     template_sim_input:
       asimmodule: geometry_optimization.atom_relax
       env_id: batch
       args:
         calculator:
           calc_id: my_mlip

----

calc\_array
-----------

**Module:** :func:`asimtools.asimmodules.workflows.calc_array.calc_array`

Run the **same asimmodule in parallel with different calculators** (or
different calculator hyperparameters). Useful for benchmarking, convergence
studies, or comparing multiple potentials.

When to use
~~~~~~~~~~~

- Compare DFT, ML potentials, and empirical force fields on the same system.
- Converge a DFT parameter (cutoff energy, k-point density) across many values.
- Benchmark a new potential against a reference.

Key parameters
~~~~~~~~~~~~~~

+-----------------------------+--------------------------------------------------+
| Parameter                   | Description                                      |
+=============================+==================================================+
| ``subsim_input``            | The sim_input for the asimmodule to run.         |
+-----------------------------+--------------------------------------------------+
| ``calculators``             | Explicit list of calculator dicts, each with     |
|                             | a ``calc_id`` or ``calc_params`` key.            |
+-----------------------------+--------------------------------------------------+
| ``template_calculator``     | Base calculator dict; ``key_sequence`` is swept  |
|                             | inside its parameters.                           |
+-----------------------------+--------------------------------------------------+
| ``key_sequence``            | Path into the calculator dict to sweep.          |
+-----------------------------+--------------------------------------------------+
| ``array_values`` /          | Values to sweep (same semantics as               |
| ``linspace_args`` /         | ``sim_array``).                                  |
| ``arange_args``             |                                                  |
+-----------------------------+--------------------------------------------------+
| ``calc_ids``                | *Deprecated* — list of calc_id strings.          |
|                             | Auto-converted to ``calculators``.               |
+-----------------------------+--------------------------------------------------+
| ``template_calc_id``        | *Deprecated* — base calc_id string.              |
|                             | Auto-converted to ``template_calculator``.       |
+-----------------------------+--------------------------------------------------+

Example — benchmark multiple calculators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   asimmodule: workflows.calc_array
   env_id: batch
   workdir: results/calc_benchmark
   args:
     calculators:
       - calc_id: emt
       - calc_id: lj_argon
       - calc_id: mace_mp
     subsim_input:
       asimmodule: singlepoint
       env_id: batch
       args:
         image: {name: Ar}
         properties: [energy, forces]

Example — converge DFT cutoff energy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   asimmodule: workflows.calc_array
   env_id: batch
   workdir: results/ecut_convergence
   args:
     template_calculator:
       calc_id: vasp_base
     key_sequence: [args, encut]
     linspace_args: [200, 600, 9]   # 200, 250, …, 600 eV
     subsim_input:
       asimmodule: singlepoint
       env_id: batch
       args:
         image: {name: Fe}
         properties: [energy]

----

distributed
-----------

**Module:** :func:`asimtools.asimmodules.workflows.distributed.distributed`

Submit an **arbitrary collection of heterogeneous sim_inputs in parallel**.
Unlike ``sim_array`` / ``image_array`` / ``calc_array``, each sub-job can be
a completely different asimmodule with different parameters.

When to use
~~~~~~~~~~~

- Run a heterogeneous set of jobs that do not share a common template.
- Fan out to many independent calculations and collect results later.
- Use as the parallel layer in a custom workflow script.

Key parameters
~~~~~~~~~~~~~~

+-----------------------------+--------------------------------------------------+
| Parameter                   | Description                                      |
+=============================+==================================================+
| ``subsim_inputs``           | Dict mapping job IDs to individual sim_inputs.   |
+-----------------------------+--------------------------------------------------+
| ``array_max``               | Max concurrent jobs in Slurm array.              |
+-----------------------------+--------------------------------------------------+
| ``group_size``              | Pack N jobs per Slurm array task.                |
+-----------------------------+--------------------------------------------------+
| ``skip_failed``             | Continue even if some sub-jobs fail.             |
+-----------------------------+--------------------------------------------------+

Example
~~~~~~~

.. code-block:: yaml

   asimmodule: workflows.distributed
   env_id: batch
   workdir: results/mixed_jobs
   args:
     subsim_inputs:
       eos_fe:
         asimmodule: eos.postprocess
         env_id: batch
         args:
           image: {name: Fe}
           calculator:
             calc_id: my_mlip
       relax_cu:
         asimmodule: geometry_optimization.cell_relax
         env_id: batch
         args:
           image: {name: Cu}
           calculator:
             calc_id: my_mlip
       sp_ar:
         asimmodule: singlepoint
         env_id: batch
         args:
           image: {name: Ar}
           calculator:
             calc_id: emt
           properties: [energy]

----

chained
-------

**Module:** :func:`asimtools.asimmodules.workflows.chained.chained`

Run asimmodules **one after the other**, where each step can depend on files
produced by the previous step. When Slurm is used, job dependencies are set
automatically via ``sbatch --dependency``.

When to use
~~~~~~~~~~~

- Multi-step workflows: relax → single-point → post-process.
- Any workflow where step N reads output files from step N-1.
- Chain a ``sim_array`` with a post-processing script.

Key parameters
~~~~~~~~~~~~~~

+-----------------------------+--------------------------------------------------+
| Parameter                   | Description                                      |
+=============================+==================================================+
| ``steps``                   | Ordered dict (``step-0``, ``step-1``, ...) of    |
|                             | sim_inputs to execute in sequence.               |
+-----------------------------+--------------------------------------------------+

.. note::

   Keys must follow the pattern ``step-N`` (zero-indexed integers). If a step
   internally launches additional Slurm jobs (e.g. via ``sim_array``), use
   ``update_dependencies`` as an intermediate step to wire up the Slurm
   dependency chain correctly.

Example — relax then compute phonons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   asimmodule: workflows.chained
   env_id: batch
   workdir: results/relax_then_phonons
   args:
     steps:
       step-0:
         asimmodule: geometry_optimization.cell_relax
         env_id: batch
         args:
           calculator:
             calc_id: my_mlip
           image: {name: Cu}
       step-1:
         asimmodule: phonons.ase_phonons
         env_id: batch
         args:
           calculator:
             calc_id: my_mlip
           image:
             image_file: ../step-0/final.xyz

Example — three-step pipeline with a parallel middle step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   asimmodule: workflows.chained
   env_id: batch
   workdir: results/pipeline
   args:
     steps:
       step-0:
         asimmodule: geometry_optimization.cell_relax
         env_id: batch
         args:
           calculator:
             calc_id: my_mlip
           image: {name: Al}
       step-1:
         asimmodule: workflows.sim_array
         env_id: batch
         args:
           key_sequence: [args, image, a]
           linspace_args: [3.9, 4.1, 5]
           template_sim_input:
             asimmodule: singlepoint
             env_id: batch
             args:
               calculator:
                 calc_id: my_mlip
               image:
                 name: Al
                 a: PLACEHOLDER
       step-2:
         asimmodule: eos.postprocess
         env_id: batch
         args:
           workdir: ../step-1

----

iterative
---------

**Module:** :func:`asimtools.asimmodules.workflows.iterative.iterative`

Run the **same asimmodule sequentially** over a sweep of values where each
job must complete before the next begins. Unlike ``sim_array``, jobs run
one at a time. Optionally, a file produced by each step can be fed as input
to the next (``dependent_file``).

When to use
~~~~~~~~~~~

- Each step's output is input to the next (e.g. active-learning loops).
- A sequential sweep where order matters.
- Iterative fitting or refinement workflows.

Key parameters
~~~~~~~~~~~~~~

+-----------------------------------+----------------------------------------------+
| Parameter                         | Description                                  |
+===================================+==============================================+
| ``template_sim_input``            | Base sim_input for each step.                |
+-----------------------------------+----------------------------------------------+
| ``key_sequence``                  | Key path swept across steps.                 |
+-----------------------------------+----------------------------------------------+
| ``array_values`` / ``linspace_args`` | Values to iterate over.                   |
| / ``arange_args``                 |                                              |
+-----------------------------------+----------------------------------------------+
| ``dependent_file``                | Filename (relative to each step's workdir)   |
|                                   | to pass as input to the next step.           |
+-----------------------------------+----------------------------------------------+
| ``dependent_file_key_sequence``   | Key path in the next step's sim_input where  |
|                                   | the file path is injected.                   |
+-----------------------------------+----------------------------------------------+

Example — sequential geometry relaxations at increasing pressures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   asimmodule: workflows.iterative
   env_id: batch
   workdir: results/pressure_series
   args:
     key_sequence: [args, pressure]
     array_values: [0, 5, 10, 20, 50]
     dependent_file: final.xyz
     dependent_file_key_sequence: [args, image, image_file]
     template_sim_input:
       asimmodule: geometry_optimization.cell_relax
       env_id: batch
       args:
         calculator:
           calc_id: my_mlip
         image: {name: Fe}
         pressure: 0

----

update\_dependencies
--------------------

**Module:**
:func:`asimtools.asimmodules.workflows.update_dependencies.update_dependencies`

An **internal helper** used by ``chained`` to wire up Slurm job dependencies
when one step internally launches additional Slurm jobs (e.g. via
``sim_array``). It reads the job IDs written by the previous step and calls
``scontrol update`` to add the correct ``afterok`` dependency on the next
step's jobs.

Most users will not need to call this directly; it is inserted automatically
by ``chained`` when needed. It has no effect outside of a Slurm environment.

.. code-block:: yaml

   # Example use inside a chained workflow (advanced)
   step-1:
     asimmodule: workflows.update_dependencies
     env_id: batch
     args:
       prev_step_dir: ../step-0
       next_step_dir: ../step-2
       skip_failed: false

----

API reference
-------------

.. toctree::
   :maxdepth: 1

   asimtools.asimmodules.workflows
