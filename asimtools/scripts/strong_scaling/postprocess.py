'''.Xauthority'''

from typing import Dict
from datetime import datetime
from glob import glob
from pathlib import Path
import numpy as np
from matplotlib import (
    pyplot as plt,
    ticker as mticker
)
from asimtools.utils import (
    write_csv_from_dict,
    read_yaml
)

def postprocess(
    outdir_pattern: str,
    outfile: str = 'output.yaml',
    job_script_file: str = 'job.sh',
    log_x: bool = True,
    log_y: bool = True,
) -> Dict:
    ''' Postprocess scaling data '''

    outdirs = glob(outdir_pattern)
    outdirs = [Path(outdir) for outdir in outdirs]
    dt_seconds = []
    ntasks = []
    nnodes = []
    for outdir in outdirs:
        outf = outdir / outfile
        job_script = outdir / job_script_file

        with open(job_script, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        for line in lines:
            if '-n ' in line:
                ntasks.append(int(line.split(' ')[-1]))
            if '-N ' in line:
                nnodes.append(int(line.split(' ')[-1]))

        output = read_yaml(outf)

        start = datetime.strptime(output['start_time'], '%Y-%m-%d %H:%M:%S.%f')
        end = datetime.strptime(output['end_time'], '%Y-%m-%d %H:%M:%S.%f')
        deltat = end - start
        dt_seconds.append(deltat.seconds)

    # Make sure ntasks are sorted from lowest to highest
    zipped = zip(ntasks, nnodes, dt_seconds)
    zipped_array = np.array(sorted(zipped, key = lambda x: x[0]))

    ntasks = zipped_array[:,0]
    nnodes = zipped_array[:,1]
    dt_seconds = zipped_array[:,2]
    speedup = dt_seconds[0] / dt_seconds

    write_csv_from_dict('scaling.csv', {
        'ntasks': ntasks,
        'nnodes': nnodes,
        'dt_seconds': dt_seconds,
        'speedup': speedup,
    })

    expected_speedup = ntasks / ntasks[0]
    fig, ax = plt.subplots()
    ax.plot(ntasks, speedup, 'k-*', label='scaling')
    ax.plot(
        ntasks,
        expected_speedup,
        'r--',
        label=f'Ideal, $\Delta t_0$(s)={dt_seconds[0]:.0f}s'
    )
    ax.set_xlabel('ntasks')
    ax.set_ylabel('Speedup')
    if log_x:
        ax.set_xscale('log', base=2)
    if log_y:
        ax.set_yscale('log', base=2)
    ax.set_xticks(ntasks)
    ax.set_xticklabels([str(ntask) for ntask in ntasks])
    ax.set_yticks(expected_speedup)
    ax.set_yticklabels([str(s) for s in expected_speedup])
    ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax.legend()

    plt.savefig('speedup_vs_tasks.png')
    plt.close(fig)

    expected_speedup = nnodes / nnodes[0]
    fig, ax = plt.subplots()
    ax.plot(nnodes, speedup, 'k-*', label='scaling')
    ax.plot(
        nnodes,
        expected_speedup,
        'r--',
        label=f'Ideal, $\Delta t_0$(s)={dt_seconds[0]:.0f}s'
    )
    ax.set_xlabel('nnodes')
    ax.set_ylabel('Speedup')
    if log_x:
        ax.set_xscale('log', base=2)
    if log_y:
        ax.set_yscale('log', base=2)
    ax.set_xticks(nnodes)
    ax.set_xticklabels([str(nnode) for nnode in nnodes])
    ax.set_yticks(expected_speedup)
    ax.set_yticklabels([str(s) for s in expected_speedup])
    ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax.legend()
    plt.savefig('speedup_vs_nodes.png')
    plt.close(fig)

    results = {
        'files': {'scaling': 'scaling.csv'}
    }

    return results
