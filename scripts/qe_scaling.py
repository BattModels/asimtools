#!/usr/bin/env python
'''
Calculate core and memory scaling of a particular system

Author: mkphuthi@github.com
'''

#pylint: disable=too-many-arguments
from typing import Sequence
import sys
from datetime import datetime
from parse_tools import get_str_btn
import numpy as np
from asimtools.job import UnitJob, Job, check_jobs
import matplotlib.pyplot as plt
from asimtools.utils import (
    write_csv_from_dict,
    get_atoms,
    join_names,
    parse_command_line,
)

def get_run_stats(outf: str, prog: str = 'PWSCF'):
    '''
    Gets the runtime, memory etc. of a calculation

    Parameters
    ----------
        outf (str): Output file path
        prog (str): Which program is being used, must be PWSCF or PHONON

    Returns
    -------
        Run time in desired format
    '''
    assert prog in ['PWSCF', 'PHONON'], 'prog must be PWSCF or PHONON'
    with open(outf, 'r', encoding='utf-8') as f:
        contents = f.readlines()

    results = {}
    for _, line in enumerate(contents):
        if line.startswith('     Number of MPI'):
            results['MPI processes'] = float(get_str_btn(
                line, ':', None
            ))
        elif line.startswith('     MPI processes distributed'):
            results['nodes'] = float(get_str_btn(
                line, 'on', 'nodes'
            ))
        elif line.startswith('     K-points division'):
            results['npool'] = float(get_str_btn(
                line, '=', None
            ))
        elif line.startswith('     Threads/MPI'):
            results['Threads/process'] = float(get_str_btn(
                line, ':', None
            ))
        elif line.startswith('     Estimated max'):
            if 'MB' in line:
                unit = 'MB'
                factor = 0.001
            else:
                unit = 'GB'
                factor = 1
            results['RAM per process'] = float(get_str_btn(
                line, '>', unit
            )) * factor
        elif line.startswith('     convergence has'):
            results['SCF iterations'] = int(get_str_btn(
                line, 'achieved in', 'iterations'
            ))
        elif line.startswith('     convergence has'):
            results['SCF iterations'] = int(get_str_btn(
                line, 'achieved in', 'iterations'
            ))
        elif line.startswith(f'     {prog}'):
            cstr = get_str_btn(line, ':', 'CPU').replace(' ', '')
            wstr = get_str_btn(line, 'CPU', 'WALL').replace(' ', '')

            for ttype, tstr in [('CPU', cstr), ('WALL', wstr)]:
                fstr = ''
                for char in tstr:
                    if char == 'd':
                        fstr += '%dd'
                    if char == 'h':
                        fstr += '%Hh'
                    if char == 'm':
                        fstr += '%Mm'
                    if char == 's' and '.' not in tstr:
                        fstr += '%Ss'
                    if char == '.':
                        fstr +='%S.%fs'

                deltat = datetime.strptime(tstr, fstr) - datetime.strptime('0', '%S')
                results[ttype] = deltat.total_seconds()

    return results

def qe_scaling(
    calc_input: dict,
    image: dict,
    nprocs_list: Sequence,
    prefix: str = '',
    workdir: str = '.',
    **kwargs
) -> dict:
    '''
    Get runtimes and memory usage
    '''
    job = Job(calc_input, {'prefix': prefix, 'workdir': workdir})
    job.start()

    # Fail early principle. There should usually be assert statements here
    # to make sure inputs are correct

    atoms = get_atoms(**image)

    sp_jobs = []
    for _, nprocs in enumerate(nprocs_list):
        sp_calc_input = calc_input.copy()
        sp_calc_input['slurm']['flags'].append(f'-n {nprocs}')
        sp_calc_input['slurm']['flags'].append(f'-J scale{nprocs}')
        nprefix = join_names([prefix, f'n{nprocs}'])
        sp_sim_input = kwargs
        atoms_input_file = join_names([nprefix, 'input_atoms.xyz'])
        sp_sim_input.update({
            'script': 'singlepoint.py',
            'image': {'input_file': atoms_input_file},
            'prefix': nprefix,
        })

        unitjob = UnitJob(
            sp_calc_input,
            sp_sim_input,
            nprefix,
        )
        unitjob.gen_input_files()
        atoms.write(unitjob.get_workdir() / atoms_input_file)
        sp_jobs.append(unitjob)

    if kwargs.get('submit', True):
        for sp_job in sp_jobs:
            sp_job.submit()
    # # Figure out if all the jobs are finished
    # all_jobs_finished = check_jobs(sp_jobs)

    # # If they are not finished, print the statues of each job
    # if not all_jobs_finished:
    #     sys.exit()

    # If all jobs are finished, prepare and write of results
    results = {
        'CPU': [], 'WALL': [], 'RAM per process': [], 'MPI processes': []
    }

    for sp_job in sp_jobs:
        # Extract run time from espresso output
        if sp_job.get_output().get('status', 'discard') != 'discard':
            stats = get_run_stats(sp_job.get_workdir() / 'QuantumEspresso.pwo')
            for key, result in results.items():
            # try:
                result.append(stats[key])
            # except:
            #     print(key, stats, result)



    table_file = join_names([prefix, 'qe_scaling_output.csv'])
    data = write_csv_from_dict(table_file, results)

    job.add_output_files({'qe_scaling': table_file})

    nprocs = np.array(data['MPI processes'])
    _, ax = plt.subplots()
    wall_speedup = np.array(data['WALL'][0]/data['WALL'])
    cpu_speedup = np.array(data['CPU'][0]/data['CPU'])
    print(nprocs, cpu_speedup)
    ax.loglog(nprocs, cpu_speedup, 'ro', label='CPU')
    ax.loglog(nprocs, wall_speedup, 'b*', label='WALL')
    ax.loglog(
        nprocs,
        np.array(nprocs),
        label='Ideal'
    )
    ax.set_xlabel('Num. MPI processes')
    ax.set_ylabel('Speedup')
    ax.set_xticks(nprocs)
    ax.set_yticks(nprocs)
    ax.set_xticklabels([str(nproc) for nproc in nprocs], rotation=90)
    ax.set_yticklabels([str(nproc) for nproc in nprocs])
    ax.set_title(f'Serial wall time = {data["WALL"][0]:.0f}s')
    ax.legend()
    plt.minorticks_off()
    plt.savefig(join_names([prefix, 'strong_scaling.png']), bbox_inches='tight')

    _, ax = plt.subplots()
    ax.loglog(nprocs, np.array(data['RAM per process']), 'o')
    ax.set_xlabel('Num. of MPI Processes')
    ax.set_ylabel('Estimated memory per process (GB)')
    ax.set_xticks(nprocs)
    ax.set_yticks(data['RAM per process'].values[0]/nprocs)
    ax.set_yticklabels([f'{val:.2f}' for val in data['RAM per process'][0]/nprocs])
    ax.set_xticklabels([str(nproc) for nproc in nprocs], rotation=90)
    plt.minorticks_off()
    plt.savefig(join_names([prefix, 'memory_scaling.png']), bbox_inches='tight')

    job.complete()

    results.update(job.get_output())
    return results


def main(argv):
    ''' Main '''
    calc_input, sim_input = parse_command_line(argv)
    qe_scaling(
        calc_input,
        **sim_input,
    )


if __name__ == "__main__":
    main(sys.argv[1:])
