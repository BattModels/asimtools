#!/usr/bin/env python
'''
Recursively check progress of jobs starting from a specified sim_input.yaml
and pring out the progress
'''
import sys
import os
from pathlib import Path
import argparse
from typing import Dict, Tuple
from colorama import Fore
from asimtools.utils import read_yaml
from asimtools.job import load_job_from_directory


def parse_command_line(args) -> Tuple[Dict, os.PathLike]:
    ''' Parse command line input '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'sim_input_file',
        metavar='sim_input_file',
        type=str,
        help='simulation input yaml file'
    )
    parser.add_argument(
        '-m',
        '--max_level',
        type=int,
        default=-1,
        help='Max. number of levels deep to print'
    )

    args = parser.parse_args(args)
    sim_input = read_yaml(args.sim_input_file)
    max_level = args.max_level
    rootdir = Path('/'.join(args.sim_input_file.split('/')[:-1]))
    return sim_input, rootdir, max_level

def main(args=None) -> None:
    """Main

    :param args: cmdline args, defaults to None
    :type args: _type_, optional
    """    
    sim_input, rootdir, max_level = parse_command_line(args)
    workdir = sim_input.get('workdir', 'results')
    if not workdir.startswith('/'):
        workdir = rootdir / workdir

    jobtree = load_job_tree(workdir)
    print_job_tree(jobtree, max_level=max_level)

def get_subjobs(workdir):
    """Get all the directories with jobs in them

    :param workdir: _description_
    :type workdir: _type_
    :return: _description_
    :rtype: _type_
    """
    subjob_dirs = []
    for subdir in workdir.iterdir():
        subsim_input = subdir / 'sim_input.yaml'
        if subsim_input.exists():
            subjob_dirs.append(subdir)
    return subjob_dirs

def load_job_tree(
    workdir: str = './',
) -> Dict:
    """Loads all the jobs in a directory in a tree format using recursion

    :param workdir: root directory from which to recursively find jobs, defaults to './'
    :type workdir: str, optional
    :return: dictionary mimicking the job tree
    :rtype: Dict
    """
    workdir = Path(workdir).resolve()

    subjob_dirs = get_subjobs(workdir)
    if len(subjob_dirs) > 0:
        subjob_dict = {}
        for subjob_dir in subjob_dirs:
            subjob_dict[subjob_dir.name] = load_job_tree(subjob_dir)
    else:
        subjob_dict = None

    job = load_job_from_directory(workdir)
    job_dict = {
        'workdir_name': workdir.name,
        'job': job,
        'subjobs': subjob_dict,
    }
    return job_dict

def get_status_and_color(job):
    ''' Helper to get printing colors '''
    status = job.get_status()[1]
    if status == 'complete':
        color = Fore.GREEN
    elif status == 'failed':
        color = Fore.RED
    elif status == 'started':
        color = Fore.BLUE
    else:
        color = Fore.WHITE
    return status, color

def print_job_tree(
    job_tree: Dict,
    indent_str: str = '',
    level: int = 0,
    max_level: int = -1,
) -> None:
    """Prints the job tree given the job_tree dictionary

    :param job_tree: Job tree to print
    :type job_tree: Dict
    :param indent_str: str to prepend, defaults to ''
    :type indent_str: str, optional
    :param level: Current printing level, defaults to 0
    :type level: int, optional
    :param max_level: Max level to print
    :type max_level: int
    """
    reset = Fore.WHITE
    subjobs = job_tree['subjobs']

    if (max_level != -1) and (level > max_level):
        pass
    elif subjobs is not None:
        workdir = job_tree['workdir_name']
        status, color = get_status_and_color(job_tree['job'])
        asimmodule = job_tree['job'].sim_input['asimmodule']
        print(color + f'{indent_str}{workdir}, asimmodule: {asimmodule},' + \
            f'status: {status}' + reset)
        if level > 0:
            indent_str = '| ' + ' ' * level
        for subjob_id in subjobs:
            print_job_tree(
                subjobs[subjob_id],
                indent_str=indent_str+'|-',
                level=level+1,
                max_level=max_level,
            )
            subjob = job_tree['job']
    else:
        subjob_dir = job_tree['workdir_name']
        subjob = job_tree['job']
        asimmodule = subjob.sim_input['asimmodule']
        status, color = get_status_and_color(subjob)
        print(color + f'{indent_str}{subjob_dir}, asimmodule: {asimmodule}, '+\
        f'status: {status}' + reset)

if __name__ == "__main__":
    main(sys.argv[1:])
