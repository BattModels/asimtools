Using built-in workflow tools
=============================

ASIMTools offers some standard tools for performing common workflows. These
are: 

#. :func:`asimtools.asimmodules.sim_array.sim_array` - Run the same asimmodule
   with one or more specified arguments of the asimmodule iterated over. This
   is the most useful asimmodule and can substitute most others except
   ``chained``

#. :func:`asimtools.asimmodules.image_array.image_array` - Run the same
   asimmodule on multiple images, e.g. a repeat calculation on a database

#. :func:`asimtools.asimmodules.calc_array.calc_array` - Run the same
   asimmodule using different calc_ids or calculator parameters based on a
   template. e.g. to converge cutoffs in DFT or benchmark many force fields

#. :func:`asimtools.asimmodules.distributed.distributed` - Run multiple
   sim_inputs in parallel

#. :func:`asimtools.asimmodules.chained.chained` - Run asimmodules one after
   the other, e.g. if step 2 results depend on step 1 etc. This allows building
   multi-step workflows.

Examples for each type of workflow are given in the examples directory and
documentation can be found in :mod:`asimtools.asimmodules`. They also serve as
templates for you to build your own workflows directly using
:func:`asimtools.job.Job` objects as an advanced user.
