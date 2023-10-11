Using built-in workflow tools
=============================

ASIMTools offers some standard tools for performing common workflows. These
are: 

#. :func:`asimtools.scripts.sim_array.sim_array` - Run the same script with one
   specified argument of the script iterated over. This is the most useful
   script and captures all the others except ``chained``

#. :func:`asimtools.scripts.image_array.image_array` - Run the same script
   on multiple structures 

#. :func:`asimtools.scripts.calc_array.calc_array` - Run the same script
   using different calculators 

#. :func:`asimtools.scripts.distributed.distributed` - Run multiple
   sim_inputs in parallel 

#. :func:`asimtools.scripts.chained.chained` - Run scripts one after the
   other, e.g. if step 2 results depend on step 1 etc.

Examples for each type of workflow are given in the examples directory and
documentation can be found in :mod:`asimtools.scripts`. They also serve as templates for
you to build your own workflows directly using :func:`asimtools.job.Job`
objects as an advanced user.
