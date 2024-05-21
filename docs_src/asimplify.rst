Developing Custom Asimmodules
=============================

This section will guide you through taking your own in-house simulation and
integrating it into asimtools. The process is designed to be as
straight-forward as reasonably possible. 

As an example, we will ASIMplify the calculation of equations of state as
shown in ASE. To understand the code, visit the ASE `tutorials page
<https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html>`. Ultimately, the
ASIMplified code will work with any element and any calculator and any
environment and be fully parallelized job submission in just a few steps!

**Use toy example** To make things easier, it is best to first replace the
expensive parts of your workflow with toy examples e.g. use an empirical
calculator like EMT instead of DFT, use a smaller/simpler structure, loop
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
        # Do one more calculation at the minimum and write to database:
        atoms.cell *= (v / atoms.get_volume())**(1 / 3)
        atoms.get_potential_energy()
        db.write(atoms, bm=B)


**Wrap in a function** The workhorse of ASIMTools is the asimmodule, any code
wrapped in a function that returns a dictionary can be run within the
ASIMTools framework. The easiest thing to do would be to take the code and
copy and paste it in a function which is defined inside an asimmodule with the
same name

.. code-block:: 
  
    from ase.build import bulk
    from ase.calculators.emt import EMT
    from ase.eos import calculate_eos
    from ase.db import connect

    def ase_eos():
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

Immediately as it is, this is an asimmodule can be run in ASIMTools. You can
run it with the following sim_input.yaml

.. code-block:: yaml

    asimmodule: /path/to/ase_eos.py 
    env_id: inline
    workdir: results

then call

.. code-block:: console

    asim-execute sim_input.yaml

Only a little bit more complicated than calling ``python ase_eos.py``
This asimmodule however still depends on specific structures and a specific
calculator. Let's do the easy thing first, let's make the asimmodule work with any
calculator using a simple change.

.. code-block:: 
  
    from ase.build import bulk
    from ase.eos import calculate_eos
    from ase.db import connect
    from asimtools.calculators import load_calc

    def ase_eos(
        calc_id,
    ):
        calc = load_calc(calc_id)
        db = connect('bulk.db')
        for symb in ['Al', 'Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']:
            atoms = bulk(symb, 'fcc')
            atoms.calc = calc
            eos = calculate_eos(atoms)
            v, e, B = eos.fit()  # find minimum
            # Do one more calculation at the minimu and write to database:
            atoms.cell *= (v / atoms.get_volume())**(1 / 3)
            atoms.get_potential_energy()
            db.write(atoms, bm=B)
        
        return {}

Just like that we can now run the asimmodule with any correctly configured
calculator for all the structures! We can even now run ``calc_array`` to
iterate getting the results using different calculators.

The final change we will make is to parallelize over structures as below

.. code-block:: 
  
    from ase.build import bulk
    from ase.eos import calculate_eos
    from ase.db import connect
    from asimtools.calculators import load_calc

    def ase_eos(
        image,
        calc_id,
    ):
        calc = load_calc(calc_id)
        db = connect('bulk.db')
        atoms = get_atoms(**image)
        atoms.calc = calc
        eos = calculate_eos(atoms)
        v, e, B = eos.fit()  # find minimum
        # Do one more calculation at the minimu and write to database:
        atoms.cell *= (v / atoms.get_volume())**(1 / 3)
        atoms.get_potential_energy()
        db.write(atoms, bm=B)
        
        return {}

Easy-peasy. We now have an asimmodule that works with arbitrary environment,
arbitrary calculator and arbitrary input structure (Of course the simulation
will fail if we give a bad structure/calculator for example)

We can do some final cleanup of the asimmodule so that it sends outputs to
``output.yaml`` and logs some checkpoints. Additionally, any asimmodules added
to the repository will need clear syntax highlighting and documentation.

.. code-block:: 
  
    from typing import Dict
    import logging
    from ase.eos import calculate_eos
    from ase.db import connect
    from asimtools.calculators import load_calc
    from asimtools.utils import get_atoms

    def ase_eos(
        image: Dict,
        calc_id: str,
        db_file: 'bulk.db'
    ) -> Dict:
        calc = load_calc(calc_id)
        db = connect(db_file)
        atoms = get_atoms(**image)
        atoms.calc = calc
        eos = calculate_eos(atoms)
        v, e, B = eos.fit()  # find minimum
        logging.info('Successfully fit EOS')
        # Do one more calculation at the minimu and write to database:
        atoms.cell *= (v / atoms.get_volume())**(1 / 3)
        atoms.get_potential_energy()
        db.write(atoms, bm=B)

        results = {'v': float(v), 'e': float(e), 'B': float(B)}
        return results

To run this asimmodule on an arbitrary structure say Argon with say the
LennardJones calculator, in a slurm job we can now use the following input
files.

sim_input.yaml:

.. code-block:: yaml

    asimmodule: /path/to/ase_eos.py 
    env_id: batch
    workdir: results
    args:
        image:
            builder: bulk
            name: Ar
        calc_id: lj_Ar

calc_input.yaml:

.. code-block:: yaml

    lj_Ar: 
        name: LennardJones
        module: ase.calculators.lj
        args:
            sigma: 3.54
            epsilon: 0.00802236
    emt:
        name: EMT
        module: ase.calculators.emt
        args: {}

env_input.yaml:

.. code-block:: yaml

    batch:
        mode:
            use_slurm: true
            interactive: false
        slurm: 
            flags:
                - -n 2
            precommands:
                - source ~/.bashrc
                - conda activate asimtools
    inline:
        mode:
            use_slurm: false
            interactive: true

Going back to the original problem, we wanted to run the simulation of multiple
different elements with the EMT calculator. To achieve that in parallel, we can
nest the ``ase_eos`` asimmodule in a
:func:`asimtools.asimmodules.workflows.sim_array.sim_array` asimmodule as follows

sim_input.yaml:

.. code-block:: yaml

    asimmodule: workflows.sim_array
    workdir: results
    args:
        key_sequence: ['args', 'image', 'name']
        array_values: ['Al', 'Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']
        env_ids: 'batch'
        template_sim_input:
            asimmodule: /path/to/ase_eos.py
            args:
                calc_id: emt
                image:
                    builder: bulk
                    crystalstructure: 'fcc'

To make the asimmodule easier to access without having to use the full path, you
can set the environment variable

.. code-block:: console

    export ASIMTOOLS_ASIMMODULE_DIR=/path/to/my/asimmodule/dir/

You can then move the ``ase_eos.py`` asimmodule to
``/path/to/my/asimmodule/dir/`` i.e. the asimmodule directory. This allows you
to refer to asimmodules prepended with the asimmodule dir as below

.. code-block:: yaml

    asimmodule: workflows.sim_array
    workdir: results
    args:
        key_sequence: ['args', 'image', 'name']
        array_values: ['Al', 'Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']
        env_ids: 'batch'
        template_sim_input:
            asimmodule: ase_eos/ase_eos.py
            args:
                calc_id: emt
                image:
                    builder: bulk
                    crystalstructure: 'fcc'

The above example loops over crystals for which ASE already has FCC lattice
parameters, but what if we want to loop over the species and corresponding
lattice parameters? We can either specify a list of ``images`` dictionaries as
``array_values`` or use ``secondary_array_values``. We can also explicitly tell
ASIMTools to include the array_values in the directory names in the standard
format (e.g. ``id-0000__Al__``, ``id-0001__Ni__`` etc.).

.. code-block:: yaml

    asimmodule: workflows.sim_array
    workdir: results
    args:
        key_sequence: ['args', 'image', 'name']
        array_values: ['Al', 'Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']
        labels: values
        secondary_key_sequences: 
        - ['args', 'image', 'a']
        secondary_array_values:
        - [4.0479, 3.524, 3.6149, 3.8907, 4.0853, 3.9242, 4.0782]
        env_ids: 'batch'
        template_sim_input:
            asimmodule: ase_eos/ase_eos.py
            args:
                calc_id: emt
                image:
                    builder: bulk
                    crystalstructure: 'fcc'
