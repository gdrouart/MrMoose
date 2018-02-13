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
#sys.path.append('/Library/Python/2.7/site-packages')
import mm_utilities as ut
from models import *

# TODO make a initial guess plotting
# TODO reproduce for a nufnu plot

def MC_Chains_plot(sampler, model_struct, fit_struct, light=None, AF_cut=0, layout=None, histo=False):
    """
    Plot the walkers vs iterations for each parameters in order to check convergence.

    :param sampler: object out of emcee, contains all the chain information
    :param model_struct: model structure defined as the combination of data/models
    :param fit_struct: fit parameters, variable, filenames
    :param light: option to plot lighter version of this plot with a defined amount of walkers
    :param AF_cut: option to filter the stick walkers from the sampler
    :param layout: option to present the plot with different line thickness, font size, etc
    :param histo: option to generate a histogram of the acceptance fraction of each chain of walkers
    """
    import matplotlib.gridspec as gridspec
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    # restore all defaults, and apply new style
    mpl.rcdefaults()
    # apply new style
    params = ut.format_plot_output(layout)
    plt.rcParams.update(params)
    # extra parameter for the split plots - smaller tick labels
    plt.rcParams.update({'xtick.labelsize':'small'})
    plt.rcParams.update({'ytick.labelsize':'small'})

    print
    print 'MC-Chains plot - convergence plot for each walker for each step'

    name_tab = ut.flatten_model_keyword(model_struct, 'param')
    min_tab = ut.flatten_model_keyword(model_struct, 'min')
    max_tab = ut.flatten_model_keyword(model_struct, 'max')
    dim = len(name_tab)
    if dim <= 5:
        max_row = dim
        max_col = 1
    else:
        max_row = int(np.ceil(dim/2.))
        max_col = 2

    # create a lighter version of the plot walkers only light keyword
    # picking randomly the defined number of walkers
    if light:
        nwalkers_loop = np.random.randint(fit_struct['nwalkers'], size=light)
    else:
        nwalkers_loop = range(fit_struct['nwalkers'])

    # identify chains of stuck walkers depending on AF_cut
    if AF_cut >= 0:
        ind_walkers = np.where(sampler.acceptance_fraction > AF_cut)[0]
    else:
        AF_plot = np.mean(sampler.acceptance_fraction)*0.5
        ind_walkers = np.where(sampler.acceptance_fraction > AF_plot)[0]

    fig = plt.figure()
    gs = gridspec.GridSpec(max_row, max_col)
    i_param = 0
    for i_col in range(max_col):
        for i_row in range(max_row):
            if i_param < dim:
                for i in nwalkers_loop:
                    ax = plt.subplot(gs[i_row, i_col])
                    ax.plot(sampler.chain[i, :, i_param], c='grey', lw=0.3)  # plot all
                for i in ind_walkers:
                    ax.plot(sampler.chain[i, :, i_param], c='k', lw=0.3)  # plot the "good" walkers
                plt.setp(ax.get_xticklabels(), visible=False)
                range_param = abs(max_tab[i_param]-min_tab[i_param])
                ax.set_ylim(min_tab[i_param]-0.05*range_param, max_tab[i_param]+0.05*range_param)
                ax.set_ylabel(name_tab[i_param])
                i_param += 1
        plt.setp(ax.get_xticklabels(), visible=True)  # set the last axis visible
        ax.set_xlabel('# steps')
    av_af = np.mean(sampler.acceptance_fraction)
    perc_filt = float((fit_struct['nwalkers']-len(ind_walkers)))/fit_struct['nwalkers']*100.
    #print len(ind_walkers),fit_struct['nwalkers']
    fig.text(0.5, 0.90,'{}, AAF={:.2f}, SW={:.1f}%'.format(fit_struct['source'],av_af,perc_filt), ha='center')
    gs.tight_layout(fig, rect=[0, 0.0, 1., 0.90])
    gs.update(hspace=0.05)
    fig.savefig(fit_struct['MCChains_plot'])

    if histo == True:
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        n, bins, patches = ax1.hist(sampler.acceptance_fraction, 20, facecolor='black', alpha=0.75)
        if AF_cut <= 0:
            ax1.axvline(x=np.mean(sampler.acceptance_fraction),ls=':',c='k')
            ax1.axvline(x=np.mean(sampler.acceptance_fraction)*0.5,ls='-',c='k')
        else:
            ax1.axvline(x=AF_cut,ls='-',c='k')
        ax1.set_xlabel('Acceptance Fraction')
        ax1.set_ylabel('# of chains')
        fig1.savefig(fit_struct['AF_histo'])
    print 'Done'


def corner_plot(sampler, model_struct, fit_struct, AF_cut=0, layout=None):
    """
    Plot the triangle plot for all the parameters.

    :param sampler: object from emcee, contain all chain information
    :param model_struct: model structure with data/model combinations
    :param fit_struct: fit parameters, variable, filenames
    :param AF_cut: acceptance fraction cut
    :return:
    """
    import corner
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    # restore all defaults, and apply new style
    mpl.rcdefaults()
    # apply new style
    params = ut.format_plot_output(layout)
    plt.rcParams.update(params)

    # get parameter names, recover dimension of table and cut MC chains
    name_tab = ut.flatten_model_keyword(model_struct, 'param')
    ndim = len(ut.flatten_model_keyword(model_struct, 'param'))

    if AF_cut >= 0:
        index_accep = np.where(sampler.acceptance_fraction > AF_cut)[0]
    else:
        index_accep = np.where(sampler.acceptance_fraction > np.mean(sampler.acceptance_fraction)*0.5)[0]
    
    samples = sampler.chain[index_accep, fit_struct['nsteps_cut']:, :].reshape((-1, ndim))
    perc_tab = np.array(fit_struct['percentiles']) * 0.01  # transformation to feed corner function

    print
    print 'Corner plot - marginalized posterior pdf of each parameter'

    # calling the corner package to plot
    fig = corner.corner(samples, labels=name_tab, show_titles=True, verbose=False,
                        quantiles=perc_tab, plot_contours=True, plot_density=False,
                        levels=perc_tab[1:-1], fill_contours=False, no_fill_contours=True,
                        smooth=1)

    # add the name of the source on top of the corner plot
    fig.text(0.9, 0.9, fit_struct['source'], ha='right')
    fig.savefig(fit_struct['triangle_plot'])

    print 'Done'

    
def SED_fnu_emcee_bestfit(data_struct, filter_struct, model_struct, fit_struct, layout=None):
    # TODO second axis
    # TODO redshift treatment
    """
    Plot the results of the fitting in nu - fnu taking confidence intervals into account.

    :param sampler:
    :param data_struct: data structure with data/model combinations
    :param filter_struct: filter structure. Same arrangement than the data_struct
    :param model_struct: model structure with data/model combinations
    :param fit_struct: fit parameters, variable, filenames
    :param layout: option to present the plot with different line thickness, font size, etc
    :return:
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    print
    print 'SED plot - Best fit visualisation'

    # restore all defaults, and apply new style
    mpl.rcdefaults()
    # apply new style
    params = ut.format_plot_output(layout)
    plt.rcParams.update(params)
    # extra parameter for the split plots - larger tick labels
    plt.rcParams.update({'xtick.labelsize':'large'})
    plt.rcParams.update({'ytick.labelsize':'large'})

    cmap = mpl.cm.rainbow
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    dim_plot = 200
    min_data = min([min(x['lambda0']) for x in data_struct])*0.1
    max_data = max([max(x['lambda0']) for x in data_struct])*10
    xscale = 10 ** (np.linspace(np.log10(min_data), np.log10(max_data), dim_plot))

    # loop on the arrangement: data and combination of models are the same color
    for i_arr in range(len(data_struct)):  # loop on the arrangement
        # get the detection and upper limits from data and plot them
        mask_d = data_struct[i_arr]['det_type'] == 'd'
        mask_ul = data_struct[i_arr]['det_type'] == 'u'

        # get the index of the models from first line of arrangements
        tmp_index_models = map(int, str.split(data_struct[i_arr]['component_number'][0], ','))
        y_total = np.zeros(dim_plot)

        for i_mod in tmp_index_models:  # loop on models in the arrangement
            tmp_colors = cmap(i_mod / float(len(data_struct)))

            # plot best fit
            if fit_struct['redshift'][i_mod] >= 0:
                y_bestfit = globals()[model_struct[i_mod]['func']](xscale, model_struct[i_mod]['bestfit'], fit_struct['redshift'][i_mod])
            else:
                y_bestfit = globals()[model_struct[i_mod]['func']](xscale, model_struct[i_mod]['bestfit'])
            ax1.plot(xscale, y_bestfit, color=tmp_colors, ls='-')
            y_total += y_bestfit

        # plot the total of components
        ax1.plot(xscale, y_total, color='k', ls='-', lw=0.6)
        # overplot the data
        ax1.errorbar(filter_struct[i_arr]['center'][mask_d],
                     data_struct[i_arr]['flux'][mask_d],
                     xerr=filter_struct[i_arr]['FWHM'][mask_d],
                     yerr=data_struct[i_arr]['flux_error'][mask_d],
                     ls='None', color='k', marker='D')
        ax1.plot(filter_struct[i_arr]['center'][mask_ul],
                 data_struct[i_arr]['flux'][mask_ul],
                 ls='None', color='k', marker='v')

    # setting plot axes, labels, etc
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel(r"F$_\nu}$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]")
    ax1.set_xlim(min([min(x['lambda0']) for x in data_struct]) * 0.1,
                 max([max(x['lambda0']) for x in data_struct]) * 10.)
    ax1.set_ylim(min([min(x['flux']) for x in data_struct]) * 0.1,
                 max([max(x['flux']) for x in data_struct]) * 10.)
    ax1.annotate(fit_struct['source'], xy=(0.9, 0.9), xycoords='axes fraction', horizontalalignment='right')

    # create the second axis as restframe frequency
    # TODO make a propoer second axis in restframe
    # if fit_struct['all_same_redshift'] is True and fit_struct['redshift'][0]>0.:
    #     ax1.annotate("%s%4.2f" % ('z=', fit_struct['redshift'][0]), xy=(0.9, 0.8),
    #                 xycoords='axes fraction', horizontalalignment='right')
    #     ax1xb = ax1.twiny()
    #     ax1xb.set_xscale("log")
    #     ax1xbticks = ax1.get_xticks()
    #     ax1xbticks_minor = ax1.get_xticks(minor=True)

    #     ax1xb.set_xticks(ax1xbticks / (1. + fit_struct['redshift'][0]))
    #     ax1xb.set_xticks(ax1xbticks_minor / (1. + fit_struct['redshift'][0]), minor=True)
    #     #ax1xb.set_xticklabels(ut.second_axis_nu_restframe(ax1xbticks, fit_struct['redshift'][0]))
    #     ax1xb.set_xlabel(r"Restframe Frequency [Hz]")
    #     ax1xb.set_xlim([i for i in ax1.get_xlim()])
    # else:
    #     print "No restframe axis plotted because components are at different redshift!"
    #     pass

    fig1.tight_layout()
    fig1.savefig(fit_struct['SED_fnu_plot'])

    print 'Done'


def split_SED_fnu_emcee_bestfit(data_struct, filter_struct, model_struct, fit_struct, layout=None):
    """
    Plot the results of the fitting in nu - fnu taking confidence intervals into account and splitting
    the SEDs per arrangement for the sake of clarity

    :param data_struct: data structure with data/model combinations
    :param filter_struct: filter structure. Same arrangement than the data_struct
    :param model_struct: model structure with data/model combinations
    :param fit_struct: fit parameters, variable, filenames
    :param layout: option to present the plot with different line thickness, font size, etc
    :return:
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    print
    print 'split SED plot - Best fit visualisation'

    # restore all defaults, and apply new style
    mpl.rcdefaults()
    # apply new style
    params = ut.format_plot_output(layout)
    plt.rcParams.update(params)
    # extra parameter for the split plots - smaller tick labels
    plt.rcParams.update({'xtick.labelsize':'x-small'})
    plt.rcParams.update({'ytick.labelsize':'x-small'})

    cmap = mpl.cm.rainbow
    nb_rows = int(np.ceil(len(data_struct)**0.5))
    nb_cols = int(np.round(len(data_struct)/len(data_struct)**0.5))
    #print "####"
    #print "Characteristics of the plot window"
    #print nb_rows, nb_cols, len(data_struct)
    #print "####"
    fig, axs = plt.subplots(nb_cols, nb_rows, sharex=True, sharey=True)
    fig.subplots_adjust(0.16, 0.11, 0.95, 0.93, 0, 0)

    dim_plot = 200
    min_data = min([min(x['lambda0']) for x in data_struct])*0.1
    max_data = max([max(x['lambda0']) for x in data_struct])*10
    xscale = 10 ** (np.linspace(np.log10(min_data), np.log10(max_data), dim_plot))

    # loop on the arrangement: data and combination of models are the same color
    for i_arr in range(len(data_struct)):  # loop on the arrangement
        # calculate the index of the subplot, trick to have origin from bottom left
        i_x = (nb_cols-1)-int(i_arr/nb_cols)
        #i_y = int(i_arr/nb_cols)
        i_y = int(i_arr % nb_cols)
        #print int(i_arr/nb_cols)
        print i_x, i_y, i_arr, data_struct[i_arr]['component'][0]
        
        # get the detection and upper limits from data and plot them
        mask_d = data_struct[i_arr]['det_type'] == 'd'
        mask_ul = data_struct[i_arr]['det_type'] == 'u'

        y_total = np.zeros(dim_plot)

        # get the index of the models from first line of arrangements
        tmp_index_models = map(int, str.split(data_struct[i_arr]['component_number'][0], ','))

        for i_mod in tmp_index_models:  # loop on models in the arrangement
            tmp_color = cmap(i_mod / float(len(data_struct)))

            # plot best fit
            if fit_struct['redshift'][i_mod] >= 0:
                y_bestfit = globals()[model_struct[i_mod]['func']]\
                    (xscale, model_struct[i_mod]['bestfit'], fit_struct['redshift'][i_mod])
            else:
                y_bestfit = globals()[model_struct[i_mod]['func']]\
                    (xscale, model_struct[i_mod]['bestfit'])
            axs[i_x, i_y].plot(xscale, y_bestfit, color=tmp_color, ls='-')
            y_total += y_bestfit

        # plot the total of components
        axs[i_x, i_y].plot(xscale, y_total, color='k', ls='-', lw=0.6)
        # overplot the data
        axs[i_x, i_y].errorbar(filter_struct[i_arr]['center'][mask_d],
                                       data_struct[i_arr]['flux'][mask_d],
                                       xerr=filter_struct[i_arr]['FWHM'][mask_d],
                                       yerr=data_struct[i_arr]['flux_error'][mask_d],
                                       ls='None', color='k', marker='D')
        axs[i_x, i_y].plot(filter_struct[i_arr]['center'][mask_ul],
                                   data_struct[i_arr]['flux'][mask_ul],
                                   ls='None', color='k', marker='v')

        # setting plot axes, labels, etc
        axs[i_x, i_y].set_xscale("log")
        axs[i_x, i_y].set_yscale("log")
        axs[i_x, i_y].set_xlim(min([min(x['lambda0']) for x in data_struct]) * 0.1,
                                       max([max(x['lambda0']) for x in data_struct]) * 10.)
        axs[i_x, i_y].set_ylim(min([min(x['flux']) for x in data_struct]) * 0.1,
                                       max([max(x['flux']) for x in data_struct]) * 10.)
        axs[i_x, i_y].annotate("{}".format(data_struct[i_arr]['component'][0]), xy=(0.9, 0.8),
                    xycoords='axes fraction', horizontalalignment='right', fontsize=8)

    # remove the empty subplots
    if len(data_struct) < nb_cols*nb_rows:
        for i_arr in range(len(data_struct), nb_cols*nb_rows):
            i_x = (nb_cols-1)-int(i_arr/nb_cols)
            i_y = int(i_arr % nb_cols)
            fig.delaxes(axs[i_x, i_y])
        
    # general title for axis
    fig.text(0.5, 0.02, "Frequency [Hz]", ha='center')
    fig.text(0.01, 0.5, r"F$_\nu}$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]", va='center', rotation='vertical')
    if fit_struct['all_same_redshift'] is True:
        fig.text(0.5, 0.95, "{}, z={}".format(fit_struct['source'],fit_struct['redshift'][0]), ha='center')
    else:
        fig.text(0.5, 0.95, "{}".format(fit_struct['source']), ha='center')
    fig.savefig(fit_struct['SED_fnu_splitplot'])

    print 'Done'

    
def SED_fnu_emcee_spaghetti(sampler, data_struct, filter_struct, model_struct, fit_struct, layout=None, AF_cut=True):
    """
    Plot the results of the fitting in nu - fnu as a spaghetti plot
    :param sampler: dictionary containing all the chains and fit exploration (emcee object)
    :param data_struct: data structure with data/model combinations
    :param filter_struct: filter structure. Same arrangement than the data_struct
    :param model_struct: model structure with data/model combinations
    :param fit_struct: fit parameters, variable, filenames
    :param layout: option to present the plot with different line thickness, font size, etc
    :param AF_cut: filter out the stuck walkers from the plot, see MC_Chains_plot
    :return:
    """
    import matplotlib as mpl
    from matplotlib.collections import LineCollection
    import matplotlib.pyplot as plt

    print
    print 'SED plot - Spaghetti visualisation'

    # restore all defaults, and apply new style
    mpl.rcdefaults()
    # apply new style
    params = ut.format_plot_output(layout)
    plt.rcParams.update(params)
    # extra parameter for the split plots - smaller tick labels
    plt.rcParams.update({'xtick.labelsize':'large'})
    plt.rcParams.update({'ytick.labelsize':'large'})

    cmap = mpl.cm.rainbow
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    dim_plot = 200
    min_data = min([min(x['lambda0']) for x in data_struct])*0.1
    max_data = max([max(x['lambda0']) for x in data_struct])*10
    xscale = 10 ** (np.linspace(np.log10(min_data), np.log10(max_data), dim_plot))

    # loop on the arrangement: data and combination of models are the same color
    for i_arr in range(len(data_struct)):  # loop on the arrangement
        # get the detection and upper limits from data and plot them
        mask_d = data_struct[i_arr]['det_type'] == 'd'
        mask_ul = data_struct[i_arr]['det_type'] == 'u'

        ax1.errorbar(data_struct[i_arr]['lambda0'][mask_d],
                     data_struct[i_arr]['flux'][mask_d],
                     data_struct[i_arr]['flux_error'][mask_d],
                     ls='None', color=cmap(i_arr / float(len(data_struct))), marker='D')
        ax1.plot(data_struct[i_arr]['lambda0'][mask_ul],
                 data_struct[i_arr]['flux'][mask_ul],
                 ls='None', color=cmap(i_arr / float(len(data_struct))), marker='v')

        # get the index of the models from first line of arrangements
        tmp_index_models = map(int, str.split(data_struct[i_arr]['component_number'][0], ','))

        for i_mod in tmp_index_models:  # loop on models in the arrangement
            # defining the parameter index for the randomn selection
            lb = sum([model_struct[i]['dim'] for i in range(i_mod)])
            ub = lb + model_struct[i_mod]['dim']
            tmp_color = cmap(i_mod / float(len(data_struct)))

            # identify chains of good walkers (not stuck) and flag the rest
            if AF_cut >= 0:
                ind_walkers = np.where(sampler.acceptance_fraction > AF_cut)[0]
            else:
                ind_walkers = np.where(sampler.acceptance_fraction > np.mean(sampler.acceptance_fraction)*0.5)[0]

            # create the line list (segs) and feed in linecollection (improvement of perf for large numbers)
            #dim_prob_full = (fit_struct['nsteps'] - fit_struct['nsteps_cut']) * fit_struct['nwalkers']
            dim_prob = (fit_struct['nsteps'] - fit_struct['nsteps_cut']) * len(ind_walkers)
            segs = np.zeros((dim_prob, dim_plot, 2))
            if fit_struct['redshift'][i_mod] >= 0:
                for cpt_s,i_steps in enumerate(range(fit_struct['nsteps_cut'], fit_struct['nsteps'])):
                    for cpt_w,i_walkers in enumerate(ind_walkers):
                        #print i_steps,i_walkers
                        y_spa = globals()[model_struct[i_mod]['func']]\
                            (xscale, sampler.chain[i_walkers, i_steps, lb:ub], fit_struct['redshift'][i_mod])
                        segs[(cpt_s*len(ind_walkers)+cpt_w), :, 1] = y_spa
            else:
                for cpt_s,i_steps in enumerate(range(fit_struct['nsteps_cut'], fit_struct['nsteps'])):
                    for cpt_w,i_walkers in enumerate(ind_walkers):
                        #print i_steps,i_walkers
                        y_spa = globals()[model_struct[i_mod]['func']]\
                            (xscale, sampler.chain[i_walkers, i_steps, lb:ub])
                        segs[(cpt_s*len(ind_walkers)+cpt_w), :, 1] = y_spa

            segs[:, :, 0] = xscale
            lines = LineCollection(segs, colors=tmp_color, lw=0.2, alpha=0.1)
            ax1.add_collection(lines)

            # overplot best fit
            if fit_struct['redshift'][i_mod] >= 0:
                y_bestfit = globals()[model_struct[i_mod]['func']]\
                    (xscale, model_struct[i_mod]['bestfit'], fit_struct['redshift'][i_mod])
            else:
                y_bestfit = globals()[model_struct[i_mod]['func']]\
                    (xscale, model_struct[i_mod]['bestfit'])
            ax1.plot(xscale, y_bestfit, color='k', ls='-')

        ax1.errorbar(filter_struct[i_arr]['center'][mask_d],
                     data_struct[i_arr]['flux'][mask_d],
                     xerr=filter_struct[i_arr]['FWHM'][mask_d],
                     yerr=data_struct[i_arr]['flux_error'][mask_d],
                     ls='None', color='k', marker='D')
        ax1.plot(filter_struct[i_arr]['center'][mask_ul],
                 data_struct[i_arr]['flux'][mask_ul],
                 ls='None', color='k', marker='v')

    # setting plot axes, labels, etc
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel(r"F$_\nu}$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]")
    ax1.set_xlim(min([min(x['lambda0']) for x in data_struct]) * 0.1,
                 max([max(x['lambda0']) for x in data_struct]) * 10.)
    ax1.set_ylim(min([min(x['flux']) for x in data_struct]) * 0.1,
                 max([max(x['flux']) for x in data_struct]) * 10.)
    ax1.annotate(fit_struct['source'], xy=(0.9, 0.9), xycoords='axes fraction', horizontalalignment='right')

    # create the second axis as restframe frequency
    # TODO make a proper second axis
    # if fit_struct['all_same_redshift'] is True and fit_struct['redshift'][0]>0.:
    #     ax1.annotate("%s%4.2f" % ('z=', fit_struct['redshift'][0]), xy=(0.9, 0.8),
    #                 xycoords='axes fraction', horizontalalignment='right')
    #     ax1xb = ax1.twiny()
    #     ax1xb.set_xscale("log")
    #     ax1xbticks = ax1.get_xticks()
    #     ax1xbticks_minor = ax1.get_xticks(minor=True)

    #     ax1xb.set_xticks(ax1xbticks / (1. + fit_struct['redshift'][0]))
    #     ax1xb.set_xticks(ax1xbticks_minor / (1. + fit_struct['redshift'][0]), minor=True)
    #     #ax1xb.set_xticklabels(ut.second_axis_nu_restframe(ax1xbticks, fit_struct['redshift'][0]))
    #     ax1xb.set_xticklabels(ax1xbticks)
    #     ax1xb.set_xlabel(r"Restframe Frequency [Hz]")
    #     ax1xb.set_xlim([i for i in ax1.get_xlim()])
    # else:
    #     print "No restframe axis plotted because components are at different redshift!"
    #     pass

    fig1.tight_layout()
    fig1.savefig(fit_struct['SED_fnu_spaplot'])
    print 'Done'


def split_SED_fnu_emcee_spaghetti(sampler, data_struct, filter_struct, model_struct, fit_struct, layout=None, AF_cut=0):
    """
    Plot the results of the fitting in nu - fnu as a spaghetti plot
    :param sampler: dictionary containing all the chains and fit exploration (emcee object)
    :param data_struct: data structure with data/model combinations
    :param filter_struct: filter structure. Same arrangement than the data_struct
    :param model_struct: model structure with data/model combinations
    :param fit_struct: fit parameters, variable, filenames
    :param layout: option to present the plot with different line thickness, font size, etc
    :param AF_cut: filter out the stuck walkers from the plot, see MC_Chains_plot
    :return:
    """
    import matplotlib as mpl
    from matplotlib.collections import LineCollection
    import matplotlib.pyplot as plt

    print
    print 'split SED plot - Spaghetti visualisation'

    # restore all defaults, and apply new style
    mpl.rcdefaults()
    # apply new style
    params = ut.format_plot_output(layout)
    plt.rcParams.update(params)
    # extra parameter for the split plots - smaller tick labels
    plt.rcParams.update({'xtick.labelsize':'x-small'})
    plt.rcParams.update({'ytick.labelsize':'x-small'})

    cmap = mpl.cm.rainbow
    nb_rows = int(np.ceil(len(data_struct)**0.5))
    nb_cols = int(np.round(len(data_struct)/len(data_struct)**0.5))
    fig, axs = plt.subplots(nb_cols, nb_rows, sharex=True, sharey=True)
    fig.subplots_adjust(0.16, 0.11, 0.95, 0.93, 0, 0)

    dim_plot = 200
    min_data = min([min(x['lambda0']) for x in data_struct])*0.1
    max_data = max([max(x['lambda0']) for x in data_struct])*10
    xscale = 10 ** (np.linspace(np.log10(min_data), np.log10(max_data), dim_plot))

    # loop on the arrangement: data and combination of models are the same color
    for i_arr in range(len(data_struct)):  # loop on the arrangement
        # calculate the index of the subplot, trick to have origin ofrom bottom left
        i_x = (nb_cols-1)-int(i_arr/nb_cols)
        i_y = int(i_arr%nb_cols)

        # get the detection and upper limits from data and plot them
        mask_d = data_struct[i_arr]['det_type'] == 'd'
        mask_ul = data_struct[i_arr]['det_type'] == 'u'

        # get the index of the models from first line of arrangements
        tmp_index_models = map(int, str.split(data_struct[i_arr]['component_number'][0], ','))

        for i_mod in tmp_index_models:  # loop on models in the arrangement
            # defining the parameter index for the randomn selection
            lb = sum([model_struct[i]['dim'] for i in range(i_mod)])
            ub = lb + model_struct[i_mod]['dim']
            tmp_color = cmap(i_mod / float(len(data_struct)))

            # identify chains of good walkers (not stuck) and flag the rest
            if AF_cut >= 0:
                ind_walkers = np.where(sampler.acceptance_fraction > AF_cut)[0]
            else:
                ind_walkers = np.where(sampler.acceptance_fraction > np.mean(sampler.acceptance_fraction)*0.5)[0]

            # create the line list (segs) and feed in linecollection (improvement of perf for large numbers)
            dim_prob = (fit_struct['nsteps'] - fit_struct['nsteps_cut']) * len(ind_walkers)
            segs = np.zeros((dim_prob, dim_plot, 2))

            if fit_struct['redshift'][i_mod] >= 0:
                for cpt_s, i_steps in enumerate(range(fit_struct['nsteps_cut'], fit_struct['nsteps'])):
                    for cpt_w, i_walkers in enumerate(ind_walkers):
                        y_spa = globals()[model_struct[i_mod]['func']]\
                            (xscale, sampler.chain[i_walkers, i_steps, lb:ub], fit_struct['redshift'][i_mod])
                        segs[(cpt_s*len(ind_walkers)+cpt_w), :, 1] = y_spa
            else:
                for cpt_s, i_steps in enumerate(range(fit_struct['nsteps_cut'], fit_struct['nsteps'])):
                    for cpt_w, i_walkers in enumerate(ind_walkers):
                        y_spa = globals()[model_struct[i_mod]['func']]\
                            (xscale, sampler.chain[i_walkers, i_steps, lb:ub])
                        segs[(cpt_s*len(ind_walkers)+cpt_w), :, 1] = y_spa
            segs[:, :, 0] = xscale
            lines = LineCollection(segs, colors=tmp_color, lw=0.2, alpha=0.1)
            axs[i_x, i_y].add_collection(lines)

            # overplot best fit
            if fit_struct['redshift'][i_mod] >= 0:
                y_bestfit = globals()[model_struct[i_mod]['func']]\
                    (xscale, model_struct[i_mod]['bestfit'], fit_struct['redshift'][i_mod])
            else:
                y_bestfit = globals()[model_struct[i_mod]['func']]\
                    (xscale, model_struct[i_mod]['bestfit'])
            axs[i_x, i_y].plot(xscale, y_bestfit, color='k', ls='-')

        # overplot the data
        axs[i_x, i_y].errorbar(filter_struct[i_arr]['center'][mask_d],
                                       data_struct[i_arr]['flux'][mask_d],
                                       xerr=filter_struct[i_arr]['FWHM'][mask_d],
                                       yerr=data_struct[i_arr]['flux_error'][mask_d],
                                       ls='None', color='k', marker='D')
        axs[i_x, i_y].plot(filter_struct[i_arr]['center'][mask_ul],
                                   data_struct[i_arr]['flux'][mask_ul],
                                   ls='None', color='k', marker='v')

        # setting plot axes, labels, etc
        axs[i_x, i_y].set_xscale("log")
        axs[i_x, i_y].set_yscale("log")
        axs[i_x, i_y].set_xlim(min([min(x['lambda0']) for x in data_struct]) * 0.1,
                                       max([max(x['lambda0']) for x in data_struct]) * 10.)
        axs[i_x, i_y].set_ylim(min([min(x['flux']) for x in data_struct]) * 0.1,
                                       max([max(x['flux']) for x in data_struct]) * 10.)
        axs[i_x, i_y].annotate("{}".format(data_struct[i_arr]['component'][0]), xy=(0.9, 0.8),
                    xycoords='axes fraction', horizontalalignment='right', fontsize=8)

    # remove the empty subplots
    if len(data_struct)<nb_cols*nb_rows:
        for i_arr in range(len(data_struct),nb_cols*nb_rows):
            i_x = (nb_cols-1)-int(i_arr/nb_cols)
            i_y = int(i_arr%nb_cols)
            fig.delaxes(axs[i_x,i_y])

    # general title for axis
    fig.text(0.5, 0.02, "Frequency [Hz]", ha='center')
    fig.text(0.01, 0.5, r"F$_\nu}$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]", va='center', rotation='vertical')
    if fit_struct['all_same_redshift'] is True:
        fig.text(0.5, 0.95, "{}, z={}".format(fit_struct['source'],fit_struct['redshift'][0]), ha='center')
    else:
        fig.text(0.5, 0.95, "{}".format(fit_struct['source']), ha='center')
    fig.savefig(fit_struct['SED_fnu_splitspaplot'])
    print 'Done'


def SED_fnu_emcee_marginalised(data_struct, filter_struct, model_struct, fit_struct, layout=None):
    """
    Plot the results of the fitting in nu - fnu taking confidence intervals into account.

    :param sampler: dictionary containing all the chains and fit exploration (emcee object)
    :param data_struct: data structure with data/model combinations
    :param filter_struct: filter structure. Same arrangement than the data_struct
    :param model_struct: model structure with data/model combinations
    :param fit_struct: fit parameters, variable, filenames
    :param layout: option to present the plot with different line thickness, font size, etc
    :return:
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    print
    print 'SED plot - Marginalised parameters visualisation'

    # restore all defaults, and apply new style
    mpl.rcdefaults()
    # apply new style
    params = ut.format_plot_output(layout)
    plt.rcParams.update(params)
    
    cmap = mpl.cm.rainbow
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    dim_plot = 200
    min_data = min([min(x['lambda0']) for x in data_struct])*0.1
    max_data = max([max(x['lambda0']) for x in data_struct])*10
    xscale = 10 ** (np.linspace(np.log10(min_data), np.log10(max_data), dim_plot))

    # generate random function for plotting in the perc intervals
    n_rand = 300

    # loop on the arrangement: data and combination of models are the same color
    for i_arr in range(len(data_struct)):  # loop on the arrangement
        # get the detection and upper limits from data and plot them
        mask_d = data_struct[i_arr]['det_type'] == 'd'
        mask_ul = data_struct[i_arr]['det_type'] == 'u'

        # get the index of the models from first line of arrangements
        tmp_index_models = map(int, str.split(data_struct[i_arr]['component_number'][0], ','))

        for i_mod in tmp_index_models:  # loop on models in the arrangement
            # defining the parameter index for the randomn selection
            #lb = sum([model_struct[i]['dim'] for i in range(i_mod)])
            #ub = lb + model_struct[i_mod]['dim']
            tmp_color = cmap(i_mod / float(len(data_struct)))

            # plot area provided by percentiles
            for i_perc in range(int(np.floor(len(model_struct[i_mod]['perc']) / 2.))):
                # generate random values in the percentile range
                rand_param = np.asarray([np.random.uniform(model_struct[i_mod]['perc'][i_perc][i_param],
                                                           model_struct[i_mod]['perc'][-(i_perc + 1)][i_param], n_rand)
                                         for i_param in range(model_struct[i_mod]['dim'])])

                list_ymod = []
                for i_rand in range(n_rand):  # loop on randomly generated models within perc
                    y_tmp = globals()[model_struct[i_mod]['func']]\
                        (xscale, rand_param[:, i_rand], fit_struct['redshift'])
                    list_ymod.append(y_tmp)
                list_ymod = np.asarray(list_ymod)
                min_list_ymod = np.min(list_ymod, axis=0)
                max_list_ymod = np.max(list_ymod, axis=0)
                ax1.fill_between(xscale, min_list_ymod, max_list_ymod,
                                 color=tmp_color, alpha=0.3 - i_perc * 0.1)

            # plot best fit
            y_bestfit = globals()[model_struct[i_mod]['func']]\
                (xscale, model_struct[i_mod]['bestfit'], fit_struct['redshift'])
            ax1.plot(xscale, y_bestfit, color=tmp_color, ls='-')

        ax1.errorbar(filter_struct[i_arr]['center'][mask_d],
                     data_struct[i_arr]['flux'][mask_d],
                     xerr=filter_struct[i_arr]['FWHM'][mask_d],
                     yerr=data_struct[i_arr]['flux_error'][mask_d],
                     ls='None', color='k', marker='D')
        ax1.plot(filter_struct[i_arr]['center'][mask_ul],
                 data_struct[i_arr]['flux'][mask_ul],
                 ls='None', color='k', marker='v')

    # setting plot axes, labels, etc
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel(r"F$_\nu}$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]")
    ax1.set_xlim(min([min(x['lambda0']) for x in data_struct]) * 0.1,
                 max([max(x['lambda0']) for x in data_struct]) * 10.)
    ax1.set_ylim(min([min(x['flux']) for x in data_struct]) * 0.1,
                 max([max(x['flux']) for x in data_struct]) * 10.)
    ax1.annotate(fit_struct['source'], xy=(0.9, 0.9), xycoords='axes fraction', horizontalalignment='right')
    ax1.annotate("%s%4.2f" % ('z=', fit_struct['redshift']), xy=(0.9, 0.8),
                 xycoords='axes fraction', horizontalalignment='right')

    # create the second axis as restframe frequency
    ax1xb = ax1.twiny()
    ax1xb.set_xscale("log")

    ax1xbticks = ax1.get_xticks()
    ax1xbticks_minor = ax1.get_xticks(minor=True)

    ax1xb.set_xticks(ax1xbticks / (1. + fit_struct['redshift']))
    ax1xb.set_xticks(ax1xbticks_minor / (1. + fit_struct['redshift']), minor=True)
    ax1xb.set_xticklabels(ut.second_axis_nu_restframe(ax1xbticks, fit_struct['redshift']))
    ax1xb.set_xlabel(r"Restframe Frequency [Hz]")
    ax1xb.set_xlim([i for i in ax1.get_xlim()])

    fig1.tight_layout()
    fig1.savefig(fit_struct['SED_fnu_margplot'])
    print 'Done'

