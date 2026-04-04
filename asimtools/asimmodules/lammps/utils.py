'''Utility functions for reading and processing LAMMPS output files.'''
from pathlib import Path
import numpy as np


def read_lammps_log(logfile, skip_failed=False):
    '''Read a LAMMPS log file and return thermodynamic data and metadata.'''
    with open(logfile, 'r', encoding='utf-8') as f:
        logtxt = f.readlines()

    starts = []
    stops = []
    natoms = None
    for i, line in enumerate(logtxt):
        if 'Step' in line:
            starts.append(i + 1)
        if len(starts) > len(stops):
            if 'Loop time' in line or 'WARNING: Pair style restartinfo' in line:
                stops.append(i)
        if 'atoms' in line:
            atoms_line = line.split()
            if len(atoms_line) == 2:
                try:
                    natoms = int(atoms_line[0])
                except ValueError:
                    pass

        if natoms is None:
            if 'Loop time' in line:
                try:
                    natoms = int(line.split()[-2])
                except ValueError:
                    pass

    if skip_failed and (len(starts) != len(stops) or len(starts) == 0):
        print(f"Incomplete run for {Path(logfile).resolve()}")
        return False

    assert len(starts) != 0, f"No data in {logfile}"
    assert natoms is not None, f"Could not find natoms in {logfile}"

    if len(starts) > len(stops):
        stops.append(-5)

    headings = logtxt[starts[0]-1].split()
    data = np.empty(len(headings))
    for start, stop in zip(starts, stops):
        for data_line in logtxt[start:stop]:
            data = np.vstack([data, np.fromstring(data_line, dtype=float, sep=' ')])

    metadata = {
        'natoms': natoms,
        'columns': headings,
    }

    return data[1:,:], metadata
