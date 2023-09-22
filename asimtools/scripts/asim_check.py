#!/usr/bin/env python
'''
Recursively check progress of jobs starting from a specified sim_input.yaml
and pring out the progress
'''
import sys
from pathlib import Path
import argparse
from typing import Dict, Tuple
from colorama import Fore
from asimtools.utils import read_yaml
from asimtools.job import load_job_from_directory


def parse_command_line(args) -> Tuple[Dict, str]:
    ''' Parse command line input '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'sim_input_file',
        metavar='simulation_input_file',
        type=str,
        help='simulation input yaml file'
    )
    args = parser.parse_args(args)
    sim_input = read_yaml(args.sim_input_file)
    return sim_input

def main(args=None) -> None:
    """Main

    :param args: cmdline args, defaults to None
    :type args: _type_, optional
    """    
    sim_input = parse_command_line(args)
    rootdir = sim_input['workdir']
    jobtree = load_job_tree(rootdir)
    print_job_tree(jobtree)

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
    workdir = Path(workdir)

    subjob_dirs = get_subjobs(workdir)
    if len(subjob_dirs) > 0:
        subjob_dict = {}
        for subjob_dir in subjob_dirs:
            subjob_dict[subjob_dir.name] = load_job_tree(subjob_dir)
    else:
        subjob_dict = None

    job = load_job_from_directory(workdir)
    job_dict = {
        'workdir': workdir.name,
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

def print_job_tree(job_tree: Dict, indent_str='', level=0) -> None:
    """Prints the job tree given the job_tree dictionary

    :param job_tree: dictionary from load_job_tree
    :type job_tree: Dict
    :param indent: indentation of output, defaults to 0
    :type indent: int, optional
    """
    reset = Fore.WHITE
    subjobs = job_tree['subjobs']
    if subjobs is not None:
        workdir = job_tree['workdir']
        status, color = get_status_and_color(job_tree['job'])
        script = job_tree['job'].sim_input['script']
        print(color + f'{indent_str}{workdir}, script: {script},' + \
            f'status: {status}' + reset)
        if level > 0:
            indent_str = '| ' + ' ' * level
        for subjob_id in subjobs:
            print_job_tree(
                subjobs[subjob_id],
                indent_str=indent_str+'|-',
                level=level+1
            )
            subjob = job_tree['job']

    else:
        subjob_dir = job_tree['workdir']
        subjob = job_tree['job']
        script = subjob.sim_input['script']
        status, color = get_status_and_color(subjob)
        print(color + f'{indent_str}{subjob_dir}, script: {script}, '+\
        f'status: {status}' + reset)

if __name__ == "__main__":
    main(sys.argv[1:])
