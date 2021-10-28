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

from . models import *
from . import mm_utilities as ut
import scipy.special
import pathos.multiprocessing as mp
from tqdm import tqdm
import emcee,ultranest
import pickle as pickle
import os.path
import logging

#create a logger
logger = logging.getLogger('mylogger')
#set logging level
logger.setLevel(logging.DEBUG)

def calc_chi2(obs_flux, error, model):
    """
    Calculated the classical chi2 likelihood value when all the data
    points are detections

    :param obs_flux: The flux given as an array containing the observed data in different bands
    :param error: Corresponding flux uncertainty, array of the same length as data
    :param model: The flux from the model at the same wavelengths
    :return: The classical chi2 value for observations without upperlimits.
    """

    return np.sum(((obs_flux-model)/error)**2)


def calc_chi2_mod(upper_lim, rms, model):
    """
    Calculates the modified chi2 likelihood from data containing upper limits
    Formula from Sawicki 2012

    :param upper_lim: Flux of the upper limit
    :param rms: The 1 sigma noise level
    :param model: Flux from the model
    :return: The modified chi2
    """
    # TODO review this formula to clean the input parameters
    return -2.*np.sum(np.log((np.pi/2.)**0.5*rms*(1.+scipy.special.erf(((rms-model)/((2**0.5)*rms))))))



def lnlike(theta, fit_struct, data, filters, param, detection_mask):
    """
    Calculated the chi2 likelihood, by calculates the model data from the initial/current guess and summing
    chi2 of each different arrangement to get the total chi2.

    :param theta: This is are the free parameters that will be fitted in the emcee process, the MCMC.
    :param fit_struct: the structure containing the information about fit parameters
    :param data: This is the observed data that is compared with the model
    :param filters: same structure than data, but to reference the filters to be used
    :param param: This contains the name of the function to be called on the different combinations of models
    :param detection_mask: List of true/false to indicate with observation is a detection or upper limit
    :return: Returns the total chi2 of the current guess.
    """
    # feeding the theta into the current parameter values
    theta = theta.flatten()
    ub = 0
    for i in range(len(param)):
        lb = ub
        ub = param[i]['dim'] + ub
        param[i]['current'] = theta[lb:ub]

    model_data = []
    number_of_component = []
    detections = []
    model_detections = []
    upper_limits = []
    model_upper_limits = []

    for i in range(len(data)):
        # make a list of list containing the number index of the model functions to use,
        # in the case of several models for
        # the same data point, the list element will be a list of numbers
        number_of_component.append(list(map(int, str.split(data[i]['component_number'][0], ','))))

        # call on the model functions via globals()
        # selects the right model (sync_lax or BB_law) via the name from parameter['func']
        # it then also sends in the data, the frequencies/wavelength and the initial or current guess
        # for the right model parameter

        min_tmp = np.log10(min(data[i]['lambda0'])*0.01)
        max_tmp = np.log10(max(data[i]['lambda0'])*100.)
        xscale = 10**np.linspace(min_tmp, max_tmp, 2000)
        
        # TODO: this is not the index i for redshift but number_of_component[i][j] to pass
        # need to figure out a way of minimising conditions and loop
        # for j
        #    if z >0
        #        global+z
        #    else
        #        global
        temp = np.zeros(2000)

        for j in range(len(number_of_component[i])):
            if fit_struct['redshift'][number_of_component[i][j]] >= 0:
                #print "pass positive, ", param[number_of_component[i][j]]['func']
                temp2 = globals()[param[number_of_component[i][j]]['func']]\
                    (xscale, param[number_of_component[i][j]]['current'], fit_struct['redshift'][number_of_component[i][j]])
            else:
                #print "pass negative, ", param[number_of_component[i][j]]['func']
                temp2 = globals()[param[number_of_component[i][j]]['func']]\
                    (xscale, param[number_of_component[i][j]]['current'])
            temp += temp2

        # Making the sum of models to go through filters
        temp_mod_filter = np.empty(data[i]['lambda0'].size)

        for j, elem in enumerate(filters[i]['name']):
            temp_mod_filter[j] = ut.integrate_filter(xscale, temp, filters[i]['wav'][j][:], filters[i]['trans'][j][:],scale=1e29)

        model_data.append(temp_mod_filter)


        # splits data in detection or upper limits, since you need to send these ones to different chi2 functions
        detections.append(data[i][detection_mask[i]])
        model_detections.append(model_data[i][np.array(detection_mask[i])])

        upper_limits.append(data[i][~detection_mask[i]])
        model_upper_limits.append(model_data[i][~np.array(detection_mask[i])])


    # calculate the total chi2 which is the main part of this function
    chi2_classic = []
    chi2_modified = []

    for i in range(len(detections)):
        chi2_classic.append(calc_chi2(detections[i]['flux'],
                                      detections[i]['flux_error'],
                                      model_detections[i]))

    for i in range(len(upper_limits)):
        chi2_modified.append(calc_chi2_mod(upper_limits[i]['flux'],
                                           upper_limits[i]['flux_error'],
                                           model_upper_limits[i]))

    return -(sum(chi2_classic)+sum(chi2_modified))


def lnprior(theta, models):
    """
    This function checks to see it the current guess it within the priors. If it is not then -inf will be returned.

    :param theta: The free parameters that are to be fitted
    :param models: Contains the min and max values of the prior
    :return: Zero or -inf, if the guess are within the allowed parameter space or not
    """

    logger.debug(theta.shape)
    # defining priors, here only uniform possible for now
    #type_func = ut.flatten_model_keyword(models, 'type')
    tmin = ut.flatten_model_keyword(models, 'min')
    tmax = ut.flatten_model_keyword(models, 'max')

    # trick to get the dimension of table right...
    theta = theta.flatten()
    # print(theta.shape)

    list_prior = [0.0 if tmin[i] < theta[i] < tmax[i] else -np.inf for i in range(len(tmin))]

    # list_prior = 0
    # for i in range(len(tmin)):
    #     if type_func[i] == 'uniform':
    #         list_prior += 0.0 if tmin[i] < theta[i] < tmax[i] else -np.inf
    #     else if type_func[i] == 'normal':
    #         list_prior += np.log10(1/(tmax[i]*(2.*np.pi)**0.5)-0.5*(param-tmin[i])**2/(tmin[i]**2)
    #     else:
    #         print "pb"
    #         sys.exit()
    
    return sum(list_prior)

def lnprior_un(theta,models):
    """
    This function check if current guess is in prior range. Ultranest implementation.
    Priors transformation in the 0 - 1 range. Only support uniform, same input file.

    :param theta: The free parameters that are to be fitted
    :param models: Contains the min and max values of the prior
    :return: Zero or -inf, if the guess are within the allowed parameter space or not
    """
    tmin = ut.flatten_model_keyword(models, 'min')
    tmax = ut.flatten_model_keyword(models, 'max')
    theta = theta.flatten()
    list_prior=np.zeros(tmin.size)
    for i,elem in enumerate(theta):
        list_prior[i] = elem*np.abs(tmax[i]-tmin[i])+tmin[i]
    return list_prior

def lnprob(theta, fit_struct, data, filters, models, detection_mask):
    # TODO implement non-uniform prior (see lnprior function)
    """
    This is the function that is called in the emcee MCMC process,it checks so that the current
    guess (theta) is within the prior and then calculated the chi2 of the current fit via lnlike.

    :param theta: The free parameters that are fitted
    :param data: Contains the observed data
    :param filters: Contains the filter used for observations
    :param models: Contains the name of which model to associate each data point
    :param detection_mask: list of true/false to decided if the data point is an detection or upperlimit
    :param redshift: redshift of the object
    :return: The maximum at posteriori estimation
    """
    lp = lnprior(theta, models)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, fit_struct, data, filters, models, detection_mask)

def decorator_like(fit_struct, data_struct, filter_struct, model_struct, detection_mask):
    """
    Decorator for the likelihood function in Ultranest
    """
    def outer(fnc):
        def inner(*args,**kwargs):
            return fnc(*args, fit_struct, data_struct, filter_struct, model_struct, detection_mask)
        return inner
    return outer

def decorator_prior(model_struct):
    """
    Decorator for the prior function in Ultranest
    """
    def outer(fnc):
        def inner(*args,**kwargs):
            return fnc(*args, model_struct)
        return inner
    return outer

def fit_source(fit_struct, data_struct, filter_struct, model_struct, Parallel=0):
    # detection mask to differenciate
    # the upper limits from detections
    # necessary to feed the chi2 calculation
    detection_mask = []
    for i in range(len(data_struct)):
        detection_mask.append(data_struct[i]['det_type'] == 'd')

    # flatten the parameters structure for processing
    ndim = len(ut.flatten_model_keyword(model_struct, 'param'))
    flat_param = ut.flatten_model_keyword(model_struct, 'current')
    flat_min = ut.flatten_model_keyword(model_struct, 'min')
    flat_max = ut.flatten_model_keyword(model_struct, 'max')
    coef_param = np.array([(flat_max - flat_min) / 2.]).flatten()

    # # initiate the walkers "ball" - TODO update with emcee function
    # if os.path.isfile(fit_struct['sampler_file']):
    #     print("file found at ",fit_struct['sampler_file'])
    #     #with open(fit_struct['sampler_file'],"rb") as f:
    #     #    chain = pickle.load(f)
    #     #last_index = chain.chain[0].nonzero()[0][-1]
    #     #print last_index
    #     #pos = [chain.chain[i][last_index] for i in range(fit_struct['nwalkers'])]
    # else:
    pos = [flat_param + 1e-4 * coef_param * np.random.randn(ndim) for i in range(fit_struct['nwalkers'])]
    print("initialise walker positions")
    last_index = 0

    print()
    print('HMC attempt for ' + fit_struct['source'])
    # single processor
#    if os.path.isfile(fit_struct['sampler_file']):
#        with open(fit_struct['sampler_file'],"rb") as f:
#            sampler = pickle.load(f)
#    else:
    # create the backend
    filename = fit_struct['sampler_file']
    backend = emcee.backends.HDFBackend(filename)
    # if restart:
    #     print('continuing last run')
    # else:
    #     backend.reset(fit_struct['nwalkers'],ndim)
    backend.reset(fit_struct['nwalkers'],ndim)

    if Parallel == 0:
        sampler = emcee.EnsembleSampler(fit_struct['nwalkers'], ndim, lnprob,
                                        args=(fit_struct, data_struct, filter_struct, model_struct, detection_mask), backend=backend)
    else:
        # multi-processing (pool created via pathos, allow to pickle the sampler)
        tmp_pool = mp.ProcessingPool(Parallel)
        sampler = emcee.EnsembleSampler(fit_struct['nwalkers'], ndim, lnprob,
                                        args=(data_struct, filter_struct, model_struct, detection_mask, fit_struct['redshift']),
                                        pool=tmp_pool, backend=backend)
            
    # progress bar (work for multiprocess or single process)
#    for sample in sampler.sample(pos,iterations=fit_source['nsteps'],progress=True):
    if fit_struct['fit_method']=='emcee': 
        with tqdm(total=fit_struct['nsteps']) as pbar:
            for i, result in enumerate(sampler.sample(pos, iterations=fit_struct['nsteps'])):
                pbar.update()
            print('HMC done!')

    #pbar = tqdm(total=fit_struct['nsteps'],initial=last_index)
    #print fit_struct['nsteps']-last_index

#    sampler.run_mcmc(None,fit_struct['nsteps'])

#    for i,_ in enumerate(sampler.sample(None, fit_struct['nsteps'])):
#        pbar.update()
#        if i % 10 == 0:
#            pass
#            with open(fit_struct['sampler_file'], 'w') as output_savefile:
#                pickle.dump(sampler, output_savefile, pickle.HIGHEST_PROTOCOL)
#                print i,' sampler saved!'
    #print('HMC done!')

    # TODO: transformation of chains
    # save the modified sampler (allows to save pools as well - pathos library allows to serialise pools)
    #with open(fit_struct['sampler_file'], 'ab') as output_savefile:
    #    pickle.dump(sampler, output_savefile, pickle.HIGHEST_PROTOCOL)
    #    print 'sampler saved!'

        return sampler
    elif fit_struct['fit_method']=='ultranest':
        param_names=list(ut.flatten_model_keyword(model_struct,'param'))

        # create the inputs functions from decorators
        @decorator_like(fit_struct, data_struct, filter_struct, model_struct, detection_mask)
        def like_func(theta,*args):
            return lnlike(theta,*args)
        @decorator_prior(model_struct)
        def prior_func(theta,*args):
            return lnprior_un(theta,*args)

        sampler_un = ultranest.ReactiveNestedSampler(param_names, like_func, prior_func)
        result=sampler_un.run()
        print('Ultranest done!')
        return result
    else:
        print('error, wrong fit_method provided')
        return 0
