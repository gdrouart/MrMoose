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

from astropy.table import Table
import os
import sys
import numpy as np
from astropy.io.votable import parse_single_table
import astropy.constants as constants


def read_data(file_name, units):
    """
    This functions read in the data file with columns
    "filter lambda0  det_type  flux   flux_error  arrangement  component"
     where the number given in "arrangement" describe which models should
     be used for each data point in the fit.

    :param units: list the units of the flux and wavelength of the input file. First element is thw flux unit,
      then wavelength unit, given a  string. e.g ['Jy','mm']
    :param file_name: The input file name
    :return: List containing the sub-tables (astropy tables), one for each
    arrangement type
    """

    # Check so that that the input file is the right format
    name, file_extension = os.path.splitext(file_name)
    flux_unit, wavelength_unit = units
    # c = 299792458  # speed of light in m/s

    # The code can only read ascii files at the moment
    if file_extension != '.dat':
        sys.exit('Input file wrong format, allowed extensions are .dat')

    if file_extension == '.dat':
        # read in the ascii table as a astropy table
        table = Table.read(file_name, format='ascii')

        # check which flux and wavelength unit have been used in the input file and
        # Convert flux to erg /s /cm2 /Hz /Sr
        # Convert wavelength to frequency in Hz
        if flux_unit == 'Jy':
            table[:]['flux'] *= 1e-23
            table[:]['flux_error'] *= 1e-23
        elif flux_unit == 'mJy':
            table[:]['flux'] *= 1e-26
            table[:]['flux_error'] *= 1e-26
        elif flux_unit == 'uJy':
            table[:]['flux'] *= 1e-29
            table[:]['flux_error'] *= 1e-29
        elif flux_unit == 'ergs':
            print 'nothing to convert, already in right units'
        else:
            sys.exit('Input flux unit not valid, allowed units are: Jy, mJy, uJy and ergs.')

        if wavelength_unit == 'm':
            table[:]['lambda0'] = constants.c.value / table[:]['lambda0']
        elif wavelength_unit == 'mm':
            table[:]['lambda0'] = np.float(constants.c.value*1e3) / table[:]['lambda0']
        elif wavelength_unit == 'um':
            table[:]['lambda0'] = np.float(constants.c.value*1e6) / table[:]['lambda0']
        elif wavelength_unit == 'Hz':
            print 'nothing to convert, already in right units'
        elif wavelength_unit == 'GHz':
            table[:]['lambda0'] *= 10**9
        elif wavelength_unit == 'MHz':
            table[:]['lambda0'] *= 10**6
        else:
            sys.exit('Input wavelength unit not valid, allowed units are: m, mm, um, Hz, MHz, GHz.')

            # Determine the number of arrangement
        number_arrangement = len(set(table['arrangement'][:]))

        # Define the table containing the subtables of the different arrangements
        sub_tables = []

        # Create a reference list containing the corresponding arrangement numbers
        refvalues = range(1, number_arrangement+1)

        for refvalue in refvalues:
            sub_tables.append(table[:][table['arrangement'] == refvalue])

    return sub_tables


def NEDtocode(filenamein, filenameout):
    """ Transform the data from NED format to readable format for MrMoose"""

    data = parse_single_table(filenamein).to_table()
    tmp = np.where((data['Frequency'] < 50e9) &
                   (data['NED_Uncertainty'] != 0.0) &
                   (data['NED_Photometry_Measurement'] != 0.0))

    # convert tables in numpy array
    # add an extra uncertainty on MSSS data
    np_passband = data['Observed_Passband'][tmp]  # name of filters
    np_passband = [np_passband[i].replace(" ", "") for i in range(np_passband.size)]
    np_freq = data['Frequency'][tmp]  # in Hz
    np_data = data['NED_Photometry_Measurement'][tmp]  # in erg/s/cm2
    for k in range(len(data['Observed_Passband'])):
        if data['Observed_Passband'][k].find('sint') == 0:
            data['NED_Uncertainty'][k] *= 10.
    np_sigma = data['NED_Uncertainty'][tmp]  # in erg/s/cm2

    np_det_type = ['d', ]*tmp[0].size
    np_comp_num = ['0', ]*tmp[0].size
    np_comp_num[-1] = '0,'
    np_comp_arr = [1, ]*tmp[0].size
    np_comments = ['"note"', ]*tmp[0].size

    table_to_write = np.array(
        zip(np_passband, np_freq, np_det_type, np_data, np_sigma, np_comp_arr, np_comments, np_comp_num),
        dtype=[('f1', 'S20'), ('f2', float), ('f3', 'S2'), ('f4', float), ('f5', float), ('f6', int), ('f7', 'S20'),
               ('f8', 'S8')])
    np.savetxt(filenameout, table_to_write, delimiter="  ",
               fmt=["%-20s"] + ["%.4e"] + ["%2s"] + ["%.4e"] * 2 + ["%2i"] + ["%-20s"] + ["%8s"],
               header="filter lambda0  det_type  flux   flux_error  arrangement  component   component_number")


def read_filters(data_struct):
    """Read the filters and store it in a structure identical to the data. """
    filter_struct = []
    for i, elem in enumerate(data_struct):
        tmp = []
        for j in elem['filter']:
            wav, trans = read_single_filter('filters/{}.fil'.format(j))
            # calculate the lambda0 and the FWHM
            center = np.average(wav, weights=trans)
            fwhm = (wav[np.where(trans > max(trans)/2.)[0][-1]] - wav[np.where(trans > max(trans)/2.)[0][0]])/2.
            tmp.append({'name': j,
                        'wav': wav,
                        'trans': trans,
                        'center': center,
                        'FWHM': fwhm})
        filter_struct.append(Table(tmp))
    return filter_struct


def read_single_filter(name):
    """ Read a single filter """
    with open(name, 'rb'):
        try:
            wav, trans = np.genfromtxt(name, unpack=True)
            return wav, trans
        except IOError:
            print "{} filter file does not exist".format(name)


def read_mod_file(name_filemod):
    """
    Read and store the models, parameters and range to use in fitting

    :param name_filemod: Name of the file containing the model details
    :return: List, of tables containing the different models and
    corresponding free parameters.
    """

    with open(name_filemod, 'rb') as f:
        list_comp = []
        func = []
        dim_func = []
        line = '#'

        # skip all the comments line
        while line.strip()[0] == '#':
            line = f.readline()
        # read the blocks of functions
        while line.strip():
            tmpfunc, dim = line.strip().split()
            func.append(tmpfunc)
            dim_func.append(dim)
            dim = int(dim)
            count = 0
            line = f.readline()
            param = []
            minparam = []
            maxparam = []
            while count < dim:
                namepar, minpar, maxpar = line.strip().split()
                param.append(namepar)
                minparam.append(minpar)
                maxparam.append(maxpar)
                count += 1
                line = f.readline()
            while line.strip() == '#':
                line = f.readline()
            minparam = np.array(map(float, minparam))
            maxparam = np.array(map(float, maxparam))
            comp = {'func': tmpfunc, 'dim': dim, 'param': param, 'min': minparam, 'max': maxparam}
            list_comp.append(comp)
    return list_comp


def set_init_guess(models):
    """
    Read in the parameters from the input model file and returns an list
    of initial guesses that are defined as the center of the provided ranges
    for the given parameters

    :param models: The list returned from function read_mod_file.
    :return: List containing tables with the initial guess for free parameter
    """

    for i in range(len(models)):
        models[i]['init'] = np.array(models[i]['min'])+(np.array(models[i]['max'])-np.array(models[i]['min']))/2.
    return models


def set_param_start(models):
    # type: (object) -> object
    """
    Set the parameter to be updated in the procedure as initial guesses

    :param models: The list returned from function read_mod_file.
    :return: List containing tables with the initial guess for free parameter
    """

    for i in range(len(models)):
        models[i]['current'] = models[i]['init']
    return models

