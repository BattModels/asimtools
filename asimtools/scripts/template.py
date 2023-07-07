#!/usr/bin/env python
'''
Describe the script briefly here. If this is the script that runs multiple steps,
describe it here. SPhinx will automatically generate docs when we get to that

Cite the papers where the method/script was first introduced here as well
Author: mkphuthi@github.com
'''

from typing import Dict # Typing is really important for our purposes
# from asimtools.calculators import load_calc #Optional: sript that uses a calculator directly needs this
from asimtools.utils import (
    # get_atoms,
    join_names,
)  #Lots of handy tools in asimtools.utils

def template(   #function name must match script name
    calc_id: str,    ## Optional: if using a calculator
    # image: Dict, # Optional: if using a structure
    # images: Dict, # Optional: if using multiple structures
) -> Dict:
    '''
    Script does xyz specifically
    '''
    
    # calc = load_calc(calc_id) # Optional: Any code that loads a calculator is replaced by this one line!
    # atoms = get_atoms(**image) # Optional
    # atoms.set_calculator(calc) # Optional

    # images = get_images(**images) # Optional

    #############################
    # Do some cool science here #
    #############################

    # We encourage a standard for naming files that describes the data and if
    # that data is for a particular condition e.g. Temperature=45K, the key value pair 
    # should appear in the name so that it's easy to know where the data is from e.g.
    # key1-value1_description_key2-value2.extension, more specifically, say you
    # run an NPT simulation for Temperature=45K then a good name is 
    # "npt_Temp-45K.thermo" with units included so that key value pairs are
    # readable, use join_names to merge parts of the name appropriately. 
    image_file = join_names([prefix, 'image_output.xyz'])
    atoms.write(image_file, format='extxyz') #extxyz is preferred format because you can store metadata

    # Collect everything you consider a "result" in a dictionary as well as 
    # the names of all files written to disk so that there is no dependence
    # on what you chose to name it.
    # Tabular data should be saved as a separate csv file and added to "files"
    # so that the results file is not clunky
    results = {
        'energy': float(energy),
        'files': {'image': image_file}
    }
    return results # You should always return a dictionary! Use {} if necessary
