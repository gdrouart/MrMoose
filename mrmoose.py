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

#import sys
#sys.path.append('/Library/Python/2.7/site-packages/')
#import numpy as np
#from astropy.cosmology import WMAP9 as cosmos
#from astropy import constants
#from guppy import hpy
import time
import cPickle as pickle
import yaml

# local package
import fitting as ft
import graphics as gp
import read_files as rd
import mm_utilities as ut
import analysis as an

from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput
from pycallgraph import Config
from pycallgraph import GlobbingFilter

# TODO double check the autocorrelation time
# TODO burn in concept good, to put back - aside for now

def SED_fit(settings_file, Parallel=None):
    """
    The main function - A to Z process for a single source and a single model file
    In case of multiple sources, see herd.py which propose a parallelisation over samples

    fit_struct contains the information about the fit
    """

    # read the files with the fit settings
    with open(settings_file, 'rb') as input:
        fit_struct = yaml.load(input)

    # some memory and timing tracking
    # h=hpy()
    # h.setref()
    # start_time=time.time()

    # expand the fit structure to include filenames, etc.
    ut.add_filenames(fit_struct)

    ### reading the data/models ###
    # input_data_file = ''.join((str.split(fit_struct['source'], '.cat5')[0], ".ascii"))
    # rd.NEDtocode(source_file, input_data_file)  # converting NED format to code format
    # data_struct = rd.readData(input_data_file, ['Jy', 'Hz'])  # reading the data

    # read data and model files
    data_struct = rd.read_data(fit_struct['source_file'], [fit_struct['unit_flux'], fit_struct['unit_obs']])
    model_struct = rd.read_mod_file(fit_struct['model_file'])

    ### reading the filter database ###
    filter_struct = rd.read_filters(data_struct)

    if fit_struct['skip_imaging'] == False:
        try:
            ut.imaging(fit_struct)
        except Exception:
            print "fits already exits, double check they are right!"
            pass
    else:
        print 'skip data imaging'
        pass

    ### init the parameter structure ###
    # Get the initial guess, which is chosen to be the median
    # value of the parameter range given in input_model_file
    # TODO get the initial parameters as two choices (on data (see theresa) or on models)
    rd.set_init_guess(model_struct)
    # Set the first parameter variables, going to be updated
    rd.set_param_start(model_struct)

    # timing
    # read_time=time.time()
    # print " time reading ---%s sec ---" % (read_time-start_time)

    ### fit the source ###
    if fit_struct['skip_fit'] == False:
        # execute fitting
        if Parallel:
            print 'multi-core sampler exploration'
            sampler = ft.fit_source(fit_struct, data_struct, filter_struct, model_struct, Parallel=Parallel)
        else:
            print 'single-core sampler exploration'
            sampler = ft.fit_source(fit_struct, data_struct, filter_struct, model_struct)
        print 'fit completed!'
    else:
        # load the sampler if not fit
        with open(fit_struct['sampler_file'], 'rb') as input:
            sampler = pickle.load(input)
        print 'sampler loaded!'

    # timing
    # fit_time=time.time()
    # print " time fitting ---%s sec ---" % (fit_time-read_time)

    # find the bestfit, percentiles, etc.
    an.find_stats_fit(sampler, model_struct, fit_struct, data_struct)

    ### plot the results ###
    layout = 'publication'
    AF_cut = 0.25   # set a value between 0 and 1, negative means a cut at mean/2
    histo=True

    # MC Chains plot to check convergence
    if fit_struct['skip_MCChains'] == False:
        gp.MC_Chains_plot(sampler, model_struct, fit_struct, layout=None, histo=histo, AF_cut=AF_cut)
    else:
        print 'skip MC Chains plotting'
        pass

    # plot the parameters confidence intervals 1D/2D
    if fit_struct['skip_triangle'] == False:
        gp.corner_plot(sampler, model_struct, fit_struct, AF_cut=AF_cut, layout=layout)
    else:
        print 'skip probability plots'
        pass

    # plot the SED with models and data
    if fit_struct['skip_SED'] == False:
        gp.SED_fnu_emcee_bestfit(data_struct, filter_struct, model_struct, fit_struct, layout=layout)
        gp.SED_fnu_emcee_spaghetti(sampler, data_struct, filter_struct, model_struct, fit_struct, layout=layout, AF_cut=AF_cut)
        #gp.SED_fnu_emcee_marginalised(data_struct, filter_struct, model_struct, fit_struct, layout=layout)

        if len(data_struct) > 1:
            gp.split_SED_fnu_emcee_bestfit(data_struct, filter_struct, model_struct, fit_struct, layout=layout)
            gp.split_SED_fnu_emcee_spaghetti(sampler, data_struct, filter_struct, model_struct, fit_struct, AF_cut=AF_cut, layout=layout)
            # gp.split_SED_fnu_emcee_marginalised(data_struct, filter_struct, model_struct, fit_struct)  #TODO
    else:
        print 'skip SED plots'
        pass

    # timing
    # plot_time=time.time()
    # print " time plotting ---%s sec ---" % (plot_time-fit_time)

    # memory use tracking
    # print h.heap()

    # calculate luminosity with uncertainties
    # save the results in a file for later use/checks
    with open(fit_struct['save_struct'], 'wb') as output:
        # format the model_struct as human readable and save
        model_sav = ut.format_sav_output(model_struct)
        yaml.dump([fit_struct,model_sav], output)

    return sampler, model_struct, data_struct, filter_struct, fit_struct

#def filtercalls(call_stack, modul, clas, func, full):
#    mod_ignore = ['scipy', 're', 'os', 'json', 'astropy', 'emcee']
#    func_ignore = []
#    clas_ignore = []
#    return modul not in mod_ignore and func not in func_ignore and clas not in clas_ignore


if __name__ == "__main__":

    graphviz = GraphvizOutput()
    graphviz.output_file = 'mrmoose_layout.png'

    config = Config()
    config.groups = False
    config.trace_filter = GlobbingFilter(exclude=[
        'pycallgraph.*',
        'tqdm.*',
        'astropy.*',
        'corner.*'
    ])

    print 'Initiating MOOSE...'
    file_example = 'fake_source_ex1.fit'
    print 'example with {}'.format(file_example)

    with PyCallGraph(output=graphviz, config=config):
        output_moose = SED_fit(file_example)

    print 'MOOSE finished!'
