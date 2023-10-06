Developing Custom Scripts
=========================

This section will guide you through taking your own in-house simulation and 
integrating it in to asimtools. The process is designed to be as 
straight-forward as reasonably possible. 

As an example, we will ASIMplify the calculation of adsorption_energies as
shown in ASE. To understand the code, visit the ASE `tutorials page <https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html>`. Ultimately, the
ASIMplified code will work with any element and , any calculator and any
environment and be fully parallelized job submission in just a few steps!

1. **Use toy example** To make things easier, it is best to first replace the
   expensive parts of your workflow with toy examples e.g. use an empirical
   calculator like EAM instead of DFT, use a smaller/simpler structure, loop
   over fewer cases etc.

The ASE example already uses a simple EMT calculator and structure which is
provided here for reference:

.. code-block:: 
    
    from ase.build import bulk
    from ase.calculators.emt import EMT
    from ase.eos import calculate_eos
    from ase.db import connect

    db = connect('bulk.db')
    for symb in ['Al', 'Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']:
        atoms = bulk(symb, 'fcc')
        atoms.calc = EMT()
        eos = calculate_eos(atoms)
        v, e, B = eos.fit()  # find minimum
        # Do one more calculation at the minimu and write to database:
        atoms.cell *= (v / atoms.get_volume())**(1 / 3)
        atoms.get_potential_energy()
        db.write(atoms, bm=B)


2. **Wrap in a function** The workhorse of ASIMTools is the function, any code
   wrapped in a function, return a dictionary can be run within the ASIMTools
   framework. The easiest thing to do would be to take the code and immediately
   put it in a function


.. code-block:: 
  
    from ase.build import bulk
    from ase.calculators.emt import EMT
    from ase.eos import calculate_eos
    from ase.db import connect

    def adsorption_energy():
        db = connect('bulk.db')
        for symb in ['Al', 'Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']:
            atoms = bulk(symb, 'fcc')
            atoms.calc = EMT()
            eos = calculate_eos(atoms)
            v, e, B = eos.fit()  # find minimum
            # Do one more calculation at the minimu and write to database:
            atoms.cell *= (v / atoms.get_volume())**(1 / 3)
            atoms.get_potential_energy()
            db.write(atoms, bm=B)
        
        return {}

.. code-block:: console

   (.venv) $ pip install lumache
