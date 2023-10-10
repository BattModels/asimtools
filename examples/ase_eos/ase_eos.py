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
