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

import mm_utilities as ut
import numpy as np
from astropy.cosmology import WMAP9 as cosmos
from astropy import constants

##############################################
# some function to use
# to be moved in the utilities
def calcul_luminosity(x, y, lim, redshift):
    # TODO make as Theresa - analysing file and move this there
    if len(lim) != 2:
        stop
    # luminosity distance
    DL = cosmos.luminosity_distance(redshift)
    dist_factor = 4. * np.pi * DL.value * DL.value * 1e4 * 1e12 * constants.pc.value * constants.pc.value

    # index corresponding to the limits of the integration -
    int_index = np.argwhere((x > (lim[0] / (1 + redshift)) * (x < (lim[1] / (1. + redshift)))))
    int_index = np.squeeze(int_index)
    y = y * (dist_factor * (1 + redshift))

    # integration
    lum = np.trapz(y[int_index], x=x[int_index])
    return lum
##############################################


def find_stats_fit(sampler, models, fit_struct, data_struct):

    # TODO here to continue updating structure of code

    ndim = len(ut.flatten_model_keyword(models, 'param'))

    # cut the chains before the provided cut
    samples = sampler.chain[:, fit_struct['nsteps_cut']:, :].reshape((-1, ndim))
    lnprob = sampler.lnprobability[:, fit_struct['nsteps_cut']:]

    # recover the index of the best fit
    bestML_index = np.array(np.unravel_index(lnprob.argmax(), lnprob.shape))
    bestML_index[1] += fit_struct['nsteps_cut']

    # extract bestfit parameters and percentiles
    best_param = sampler.chain[bestML_index[0], bestML_index[1], :]
    perc_param = np.percentile(samples, fit_struct['percentiles'], axis=0)

    # add to the model structure
    lb = 0
    for i in range(len(models)):
        ub = lb + models[i]['dim']
        models[i]['bestfit'] = best_param[lb:ub]
        models[i]['perc'] = [perc_param[j][lb:ub] for j in range(len(perc_param))]
        lb += models[i]['dim']
    # the float() is a trick to force saving the valu in a variable and not only create a pointer
    fit_struct['best_lnL'] = float(sampler.lnprobability[bestML_index[0], bestML_index[1]])

    # TODO extract the seed number for the RNG in emcee for possible future reproduction (needs emcee and init position)

    # Calculate the AICc and save
    # in developpement
    # consider each flux of each arrangement as one datapoint
    ndata = sum([len(elem['filter']) for i,elem in enumerate(data_struct)])
    try:
        print "AICc calculation... still in developpement, to be used with caution"
        penalty_factor = 2 * ndim * (ndim +1) / (ndata - ndim -1)
        fit_struct['AICc'] = 2* ndim -2 * fit_struct['best_lnL'] + penalty_factor
    except:
        print "AICc cannot be calculated, too many parameters compared to data"
