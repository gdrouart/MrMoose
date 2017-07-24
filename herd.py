"""
MrMoose: A multi-wavelength, multi-resolution, multi-sources fitting code.
Allow to simultaneously fit multiple models in order to interpret the spectral
energy distribution of an astrophysical source from simple to more complex case, 
making use of a bayesian framework.

Copyright (C) 2017 Guillaume Drouart, Theresa Falkendal

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
"""

import mrmoose as mm
import glob
import numpy as np
import pathos.multiprocessing as mp
from tqdm import tqdm

def herd(list_input):
    """ Wrapper to parallelise MrMoose on samples"""

    ndim = len(list_input)
    Parallel_cpu = mp.cpu_count()  # take the number of CPUs
    print '# files: ', ndim
    print '# cores: ', Parallel_cpu

    # create the pool of jobs
    pool = mp.ProcessingPool(Parallel_cpu)

    # run and display progress bar
    with tqdm(total=ndim):
        herd_outputs = pool.map(mm.SED_fit, list_input)
    return herd_outputs

if __name__ == "__main__":
    # looping over mr_moose
    list_file = glob.glob('data/*cat5')
    name_z = np.genfromtxt('shzrg.dat', dtype=[('name', 'S20'), ('z', 'f8')], unpack=True)

    # create and duplicate the input structure
    input_struct = {"source_file": list_file[0],
                    "redshift": name_z[0]['z'],
                    "model_file": 'inputs/models1.mod',
                    "nwalkers": 20,
                    "nsteps": 20,
                    "nsteps_cut": 18,
                    "percentiles": np.array([10., 25., 50., 75., 90.]),
                    "skip_fit": False,
                    "skip_MCChains": False,
                    "skip_triangle": False,
                    "skip_SED": False}
    ndim = len(list_file)
    list_input = [input_struct, ]*ndim

    # Change the names to the sources
    for i in range(ndim):
        list_input[i]['source_file'] = list_file[i]
        list_input[i]['redshift'] = name_z[i]['z']

    # execute the loop
    output_herd = herd(list_input)

