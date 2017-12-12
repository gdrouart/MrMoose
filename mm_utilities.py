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

import numpy as np
import itertools
from models import *

def flatten_model_keyword(models, keyword):
    """Flatten the model file for the provided keyword in 1D array """
    tmp_list=[]
    for i in range(len(models)):
        tmp_list.append(models[i][keyword])
    return np.asarray(list(itertools.chain.from_iterable(tmp_list)))


def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific notation of
    the given number formatted for use with Latex or Mathtext, with
    specified number of significant decimal digits and precision (number
    of decimal digits to show). The exponent to be used can also be
    specified explicitly.

    :param num: the number to format
    :param decimal_digits:  number of decimal digits if 0 only an integer is returned, if -1 only exponent
    :param precision: precision to reach
    :param exponent: exponent to use
    :return: a string with the right format
    """
    if not exponent:
        exponent = int(np.floor(np.log10(abs(num))))
        coeff = np.round(num / np.float(10 ** exponent), decimal_digits)
        if not precision:
            precision = decimal_digits
        if decimal_digits == -1:
            return r"$10^{{{1:d}}}$".format(coeff, exponent, precision)
        elif decimal_digits == 0:
            coeff = np.int(coeff)
            return r"${0:d}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)
        else:
            return r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)


def second_axis_nu_restframe(x, redshift):
    """ Calculate the restframe axis """
    v = x
    return [sci_notation(x, -1) for x in v]


def second_axis_nu(x):
    """ Calculate the second axis """
    v = x
    return [sci_notation(x, -1) for x in v]


def format_sav_output(mod_struct):
    """ Transform numpy array in list in the model structure to make the yaml dump human readable"""
    import copy

    mod_struct_cop = copy.deepcopy(mod_struct)
    for i in range(len(mod_struct)):
        for j in mod_struct[i].keys():
            mod_struct_cop[i][j] = np.array(mod_struct[i][j]).tolist()
    return mod_struct_cop

def format_plot_output(layout=None):
    # TODO add a precision mode - error smaller than symbols
    """ Choose a layout for the graphics """
    if layout=='publication':
        params={'font.size': 14.,
                'font.family':'serif',
                'font.weight':'bold',
                'axes.linewidth': 2.5,
                'axes.labelsize':18.,
                'axes.labelweight':'bold',
                'mathtext.fontset':'custom',
                'mathtext.bf':'serif:bold',
                'mathtext.tt':'serif',
                'xtick.major.size': 6.,
                'xtick.minor.size': 4.,
                'ytick.major.size': 6.,
                'ytick.minor.size': 4.,
                'xtick.major.width': 2.,
                'xtick.minor.width': 1.,
                'ytick.major.width': 2.,
                'ytick.minor.width': 1.,
                'lines.markersize':12,
                'lines.linewidth':2.}
    elif layout=='talk':
        params={'font.size': 18.,
                'font.family':'serif',
                'font.weight':'bold',
                'axes.linewidth': 2.5,
                'axes.labelsize':20.,
                'axes.labelweight':'bold',
                'mathtext.fontset':'custom',
                'mathtext.bf':'serif:bold',
                'mathtext.tt':'serif',
                'xtick.major.size': 6.,
                'xtick.minor.size': 4.,
                'ytick.major.size': 6.,
                'ytick.minor.size': 4.,
                'xtick.major.width': 2.,
                'xtick.minor.width': 1.,
                'ytick.major.width': 2.,
                'ytick.minor.width': 1.,
                'lines.markersize':12,
                'lines.linewidth':2.}
    else:
        params={'font.size':12.,
                'font.family':'serif'}
    return params


def add_filenames(fit_struct):
    """ Function to centralise and prepare filenames for outputs"""

    fit_struct['model'] = '.'.join(str.split(str.split(fit_struct['model_file'], '/')[-1], '.')[:-1])
    fit_struct['source'] = '.'.join(str.split(str.split(fit_struct['source_file'], '/')[-1], '.')[:-1])

    directory = 'outputs/'
    suffix = fit_struct['source'] + '_' + fit_struct['model'] \
             + '_w' + str(fit_struct['nwalkers']) + '_s' + str(fit_struct['nsteps'])

    fit_struct['triangle_plot'] = directory + suffix + '_triangle.pdf'
    fit_struct['MCChains_plot'] = directory + suffix + '_MCChains.pdf'
    fit_struct['AF_histo'] = directory + suffix + '_AF_histo.pdf'
    fit_struct['SED_fnu_plot'] = directory + suffix + '_SED_fnu.pdf'
    fit_struct['SED_fnu_splitplot'] = directory + suffix + '_SED_fnu_split.pdf'
    fit_struct['SED_fnu_spaplot'] = directory + suffix + '_SED_fnu_spag.pdf'
    fit_struct['SED_fnu_splitspaplot'] = directory + suffix + '_SED_fnu_spag_split.pdf'
    fit_struct['SED_fnu_margplot'] = directory + suffix + '_SED_fnu_marg.pdf'
    fit_struct['SED_fnu_splitmargplot'] = directory + suffix + '_SED_fnu_marg_split.pdf'  # TODO
    #fit_struct['SED_nufnu_plot'] = directory + suffix + '_SED_nufnu.pdf'
    fit_struct['SED_file'] = directory + suffix + '.sed'  # TODO
    fit_struct['sampler_file'] = directory + suffix + '.pkl'
    fit_struct['save_struct'] = directory + suffix + '.sav'
    return fit_struct


def integrate_filter(sed_nu, sed_flux, filter_nu, filter_trans):
    """
    Function integrating SED flux over a specific filter

    :param sed_nu: array of frequencies of the SED
    :param sed_flux: array of flux of the SED (at each frequency)
    :param filter_nu: array of frequencies of the filter
    :param filter_trans: array of transmission of the filter at each frequency
    :return: integrated flux
    """
    filter_trans_interp = np.interp(sed_nu, filter_nu, filter_trans, right=0.0, left=0.0)

    if min(filter_nu) < min(sed_nu) or max(filter_nu) > max(sed_nu):
        print "filter outside range!"
        return 0.0
    else:
        # find max filter
        i1 = np.where(filter_trans_interp >= max(filter_trans_interp)/5.)[0][0]
        i2 = np.where(filter_trans_interp >= max(filter_trans_interp)/5.)[0][-1]

        ic = np.floor((i1 + i2) / 2.)
        scale_w = np.abs(sed_flux[int(ic)])

        k1 = 0
        if np.isfinite(scale_w) and scale_w > 0:
            k1 = scale_w * np.trapz((sed_flux * sed_nu / scale_w) * filter_trans_interp, x=sed_nu)
        k2 = np.trapz(sed_nu * filter_trans_interp, x=sed_nu)
        return k1 / k2


def create_filter_gate(name, freq, trans_window, size=500, trans_value=1.0):
    """
    Create a fake filter of the gate form. transimission of 1 in a certain frequency range, 0 elsewhere.

    :param name: name of the filter
    :param freq: array of frequency
    :param trans_window: range of frequency to consider at the given transmission value
    :param size: dimension of the array
    :param trans_value: transmission value to assume
    """
    wav = np.linspace(min(freq), max(freq), size)
    trans = np.zeros(size)
    indexmin = np.argmin(np.abs(wav-min(trans_window)))
    indexmax = np.argmin(np.abs(wav-max(trans_window)))
    trans[indexmin:indexmax] = trans_value
    with open(name, 'wb'):
        np.savetxt(name, np.vstack((wav, trans)).T, header='Freq [Hz]   Transmission')
    return


def twoD_Gaussian((x,y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    """
    Create a 2D gaussian

    :param amplitude: maximum amplitude of the Gaussian
    :param xo: x-axis pixel coordinate
    :param yo: y-axis pixel coordinate
    :param sigma_x: standard deviation on x-axis
    :param sigma_y: standard deviation on y-axis
    :param theta: angle of the Gaussian
    :param offset: level of "background" signal
    :return: 2D array of values
    """
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g


def create_image(namefits, RA_source, Dec_source, Flux_source, res, pixsize):
    """
    Create a fake fits image corresponding to a single filter from the datafile

    :param namefits: name of the image
    :param RA_source: right ascension of the source(s)
    :param Dec_source: declination of the source(s)
    :param Flux_source: flux of the source(s)
    :param res: resolution at this frequency of observation (in arcsec)
    :param pixsize: size of the pixel (in arcsec)
    """
    from astropy import wcs
    from astropy.io import fits
    from astropy.coordinates import Angle

    # transform the coordinates, resolution and pixsize, in degrees
    RA_source = Angle(RA_source).degree
    Dec_source = Angle(Dec_source).degree
    res = [i/3600. for i in res]
    pixsize = [i/3600. for i in pixsize]

    # definition of the size of the image
    av_center_RA = np.mean(RA_source)
    av_center_Dec = np.mean(Dec_source)
    delta_RA = np.abs(min(RA_source) - max(RA_source)) + max(res) * 3.
    delta_Dec = np.abs(min(Dec_source) - max(Dec_source)) + max(res) * 3.
    dim = np.floor((max([delta_RA, delta_Dec]) + max(res)) / min(pixsize))

    # create the WCS used to model the flux
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [dim / 2,
                   dim / 2]  # reference pixel
    w.wcs.cdelt = np.array([-min(pixsize),
                            min(pixsize)])  # pixel resolution -- only min instance taken
    w.wcs.crval = [av_center_RA,
                   av_center_Dec]  # coordinate reference pixel
    w.wcs.ctype = ["RA-TAN",
                   "DEC-TAN"]  # projection

    pix_source = w.wcs_world2pix(zip(RA_source, Dec_source), 1)

    # Create x and y indices -- need to start one because no 0th pixel on image
    x = np.linspace(1, dim, dim)
    y = np.linspace(1, dim, dim)
    x, y = np.meshgrid(x, y)

    # create data
    data = sum(
        [twoD_Gaussian((x, y), Flux_source[i], pix_source[i][0], pix_source[i][1],
                       res[i] / pixsize[i] / 2, res[i] / pixsize[i] / 2, 0, 0)
         for i, elem in enumerate(Flux_source)])

    # save the header, the data in the fits file
    header = w.to_header()
    hdu = fits.PrimaryHDU(header=header, data=data)
    hdu.writeto(namefits)


def imaging(fit_struct):
    """Create the images for the entire datafile"""
    from astropy.table import Table, join
    from itertools import groupby

    # read again the table (simpler format) and sort it
    table = Table.read(fit_struct['source_file'], format='ascii')
    table.sort('filter')

    # find the same filter name
    for group, elem, in groupby(table, lambda x: x['filter']):
        tmp_ra = []
        tmp_dec = []
        tmp_flux = []
        tmp_res = []
        tmp_pixsize = []
        for filt in elem:
            tmp_ra.append(filt['RA'])
            tmp_dec.append(filt['Dec'])
            tmp_flux.append(filt['flux'] if filt['det_type'] == 'd' else 0.)
            tmp_res.append(filt['resolution'])
            tmp_pixsize.append(filt['resolution']/5.)
        tmp_fitsname = 'outputs/' + fit_struct['source'] + '_' + filt['filter'] + '.fits'
        create_image(tmp_fitsname, tmp_ra, tmp_dec, tmp_flux, tmp_res, tmp_pixsize)

    return 0


def save_bestfit_SED(data_struct, fit_struct, model_struct):
    #TODO create a file with the SED as frequency vs flux, one column for each model
    #TODO create a header to explain each column

    print "SED saving in developement, use with caution"

    tab = []
    min_data = min([min(x['lambda0']) for x in data_struct])*0.1
    max_data = max([max(x['lambda0']) for x in data_struct])*10
    xscale = 10 ** (np.linspace(np.log10(min_data), np.log10(max_data), 200))
    tab.append(xscale)
    for i_arr in range(len(data_struct)):
        tmp_index_models = map(int, str.split(data_struct[i_arr]['component_number'][0], ','))
        for i_mod in tmp_index_models:
    # overplot best fit
            if fit_struct['redshift'][i_mod] >= 0:
                y_bestfit = globals()[model_struct[i_mod]['func']]\
                    (xscale, model_struct[i_mod]['bestfit'], fit_struct['redshift'][i_mod])
            else:
                y_bestfit = globals()[model_struct[i_mod]['func']]\
                    (xscale, model_struct[i_mod]['bestfit'])
            tab.append(y_bestfit)
    # open and save the SED in file .sed
    np.savetxt(fit_struct['SED_file'],np.transpose(tab), delimiter=',', header="Freq[Hz], "+model_struct[0]['func']+"[erg/s/cm/Hz]")
    return

