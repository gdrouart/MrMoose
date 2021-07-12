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

"""
Library of models for the code to fit
Functions can be added here by the user
Each model can come in two flavours, with or without redshift constraint
This can be particularly useful when fitting objects at different redshifts
blended at some wavelengths
"""

import astropy
import numpy as np
from astropy import analytic_functions

bbfunction = analytic_functions.blackbody.blackbody_nu

# TODO have a look on optimising BB_law and AGN_law or the upper limit

def sync_law(x, param, redshift):
    # TODO beware of redshift use here
    """
    Flux_nu=norm*freq^alpha
    This is a function that describe the synchrotron radiation as a power
    law. Input is the observed frequency, which is the shifted to that the
    power law is calculated at rest frame frequency. The flux is then
    shifted to the observed fame, by taking into account the decrease of
     flux due to distance.
    :param x: Frequency in Hz
    :param param: List of the normalization factor and the power alpha, param=[norm,aplha]
    :param redshift: The redshift of the object
    :return: The observed flux of synchrotron radiation as described by a power law
    """
    norm, alpha = param
    # moving model in the observed frame
    x_z = x / (1. + redshift)
    y_z = (x_z**alpha)*(10**norm)*(1. + redshift)
    return y_z


def sync_law_z(x, param):
    """
    Flux_nu=norm*freq^alpha
    This is a function that describe the synchrotron radiation as a power
    law. Input is the observed frequency, which is the shifted to that the
    power law is calculated at rest frame frequency. The flux is then
    shifted to the observed fame, by taking into account the decrease of
    flux due to distance. This function is similar to sync_law but allow for
    the redshift as a free parameter
    :param x: Frequency in Hz
    :param param: List of the normalization factor and the power alpha and redshift, param=[norm, alpha, z]
    :return: The observed flux of synchrotron radiation as described by a power law
    """
    # function to be fully usable in version --- ALPHA #
    # TODO - to use to have redshift as free parameter, see upcoming update
    norm, alpha, z = param
    x_z = x / (1. + z)
    y_z = (x_z**alpha)*(10**norm)*(1. + z)
    return y_z


def BB_law(x, param, redshift):
    """
    Calculated the blackbody radiation from head dust. Inout is the observed
    frequency, normalizing factor and the temperature. It shifts to restframe
    frequency and calculate the blackbody spectrum. The output is the specific
    intensity which is constant over distances.
    :param x: Frequency in Hz
    :param param: List of the normalization factor and temperature, param=[norm,temp]
    :param redshift: The redshift z of the object
    :return: the observed flux_nu of the blackbody model
    """
    from astropy import analytic_functions
    bbfunction = analytic_functions.blackbody.blackbody_nu
    norm, temp = param
    redshift_shifted_freq = (1 + redshift)*x
    stradian = 4*np.pi
    
    # The astropy.blackbody_nu gives the specific intensity of the blackbody ragiation in
    # units erg s^-1 cm^-2 Hz^-1 sr^-1, the specific intenisty is
    # constant over distance due to the r^2 dependence of the solid angle and flux cancel out each other.
    # but to have a comparable quanitiy to the oberved flux we need to multiply by
    # the solid angle, which for a sphere is 4pi
    # return stradian*(10**norm)*analytic_functions.blackbody.blackbody_nu(redshift_shifted_freq, temp).value
    y_z = stradian*(10**norm)*bbfunction(redshift_shifted_freq, temp).value
    y_z *= (1. + redshift)
    return y_z


def BB_law_z(x, param):
    """
    Calculated the blackbody radiation from head dust. Inout is the observed
    frequency, normalizing factor and the temperature. It shifts to restframe
    frequency and calculate the blackbody spectrum. The output is the specific
    intensity which is constant over distances.
    :param x: Frequency in Hz
    :param param: List of the normalization factor and temperature, param=[norm,temp]
    :param redshift: The redshift z of the object
    :return: the observed flux_nu of the blackbody model
    """
    from astropy import analytic_functions
    bbfunction = analytic_functions.blackbody.blackbody_nu
    norm, temp, redshift = param
    redshift_shifted_freq = (1 + redshift) * x
    stradian = 4*np.pi
    
    # The astropy.blackbody_nu gives the specific intensity of the blackbody ragiation in
    # units erg s^-1 cm^-2 Hz^-1 sr^-1, the specific intenisty is
    # constant over distance due to the r^2 dependence of the solid angle and flux cancel out each other.
    # but to have a comparable quanitiy to the oberved flux we need to multiply by
    # the solid angle, which for a sphere is 4pi
    # return stradian*(10**norm)*analytic_functions.blackbody.blackbody_nu(redshift_shifted_freq, temp).value
    y_z = stradian*(10**norm)*bbfunction(redshift_shifted_freq, temp).value
    y_z /= (1. + redshift)
    return y_z


def MBB_law(nu_obs, param, redshift):
    """
    Calculate the flux for a modified black body. 
    Formula taken from Gilli et al. 2014

    param nu: array of frequency
    param param: array of the parameter values
    param redshift: redshift of the source
    return: calculated flux for the frequencies provided in the array with the given parameter

    parameter array description:
    p[0]: norm => normalisation fator in log
    p[1]: temp => temerature of the BB in Kelvin
    p[2]: beta => index of the slope
    p[3]: nu_0 => pivotal frequency for the two different regime of emission in log
    """
    from astropy import analytic_functions
    bbfunction = analytic_functions.blackbody.blackbody_nu
    norm, beta, temp = param
#    beta = 2.0
    nu_0 = np.log10(1.5e12)

    nu_rest = nu_obs * (1+redshift) # shift the frequency range to restframe
    y_z = (10**norm)*bbfunction(nu_rest, temp).value * (1. - np.exp(-1.0*(nu_rest/10**nu_0)**beta))
    y_z *= (1. + redshift) # shift the flux to observed frame
    return y_z

def MBB_law_z(nu_obs, param):
    """
    Calculate the flux for a modified black body. 
    Formula taken from Gilli et al. 2014

    param nu: array of frequency
    param param: array of the parameter values
    param redshift: redshift of the source
    return: calculated flux for the frequencies provided in the array with the given parameter

    parameter array description:
    p[0]: norm => normalisation fator in log
    p[1]: temp => temerature of the BB in Kelvin
    p[2]: beta => index of the slope
    p[3]: nu_0 => pivotal frequency for the two different regime of emission in log
    p[4]: redshift  => redshift
    """
    from astropy import analytic_functions
    bbfunction = analytic_functions.blackbody.blackbody_nu
    norm, temp, beta, redshift = param
#    beta = 2.0
    nu_0 = np.log10(1.5e12)

    nu_rest = nu_obs * (1+redshift) # shift the frequency range to restframe
    y_z = (10**norm)*bbfunction(nu_rest, temp).value * (1. - np.exp(-1.0*(nu_rest/10**nu_0)**beta))
    y_z *= (1. + redshift) # shift the flux to observed frame
    return y_z

def double_sync_law(nu, param, redshift):
    """
    Calculate the flux for a double power law with break frequency.
    Inspired from the SSRQ formula in Polletta et al. 2000.

    :param nu: array of increasing frequency (in log)
    :param param: array of the parameter values
    :param redshift: redshift of the source
    :return: the calculated flux at the provided frequency array with the given parameters

    parameters array description:
    p[0]: norm => normalisation factor in log
    p[1]: nu_break => break frequency between spectral index in log
    p[2]: alpha1 => lower frequency index
    p[3]: alpha2 => higher frequency index
    """
    n, nu_break, a1, a2 = param
    y = np.float128(10**n*(nu*(1.+redshift)/10**nu_break)**a1) * \
        (1.-np.float128(np.exp(-1.*(10**nu_break/(nu*(1.+redshift)))**(a1-a2))))
    y /= (1+redshift)
#    y /= nu teporary hack
    return np.float64(y)

def double_power_law(nu,a,redshift,s=1.):
    norm,nu_t,a1,a2=a
    y = norm*nu**a1*(1.+(nu/(10**nu_t))**(np.abs(a1-a2)*s))**(np.sign(a2-a1)/s)
    y /= (1+redshift)
    return y

def double_sync_law_z(nu, param):
    """
    Calculate the flux for a double power law with break frequency.
    Inspired from the SSRQ formula in Polletta et al. 2000.

    :param nu: array of increasing frequency (in log)
    :param param: array of the parameter values
    :return: the calculated flux at the provided frequency array with the given parameters

    parameters array description:
    p[0]: norm => normalisation factor in log
    p[1]: nu_break => break frequency between spectral index in log
    p[2]: alpha1 => lower frequency index
    p[3]: alpha2 => higher frequency index
    p[4]: redshift 
    """
    n, nu_break, a1, a2, z = param
    nu_z = nu / (1 + z)
    y = np.float128(10**n*(nu*(1.+z)/10**nu_break)**a1) * \
        (1.-np.float128(np.exp(-1.*(10**nu_break/(nu*(1.+z)))**(a1-a2))))
    y /= (1+z)
    y /= nu
    return np.float64(y)


def absorp_double_sync_law(nu, param, redshift):
    """
    Calculate the flux for a double power law and absorption at low frequency.
    Inspired from Polletta et al. 200 formula. Contains 3 spectral slopes and 2 break frequencies

    :param nu: array of increasing frequency (in log)
    :param param: array of the parameter values
    :param redshift: redshift of the source
    :return: the calculated flux at the provided frequency array with the given parameters

    parameters array description:
    p[0]: norm => normalisation factor in log
    p[1]: nu_to => frequency where tau=1 (optically thick/thin transition) in log
    p[2]: nu_break => break in higher frequency spectral index in log
    p[3]: alpha0 => absorption index
    p[4]: alpha1 => lower frequency index
    p[5]: alpha2 => higher frequency index
    """
    n,nu_to,nu_break,a0,a1,a2=param
#    term1 = np.float128((nu*(1.+redshift)/(10**nu_break))**-a1)
#    term2 = 1.-np.float128(np.exp(-1.*(10**nu_break/(nu*(1.+redshift)))**(a1-a2)))
#    term3 = 1.-np.float128(np.exp(-1.*(10**nu_to/(nu*(1.+redshift)))**(a0-a1)))
#    y = 10**n*term1*term2#/term3
    y = 10**n * np.float128((nu/10**nu_break)**a1) * \
        (1.-np.float128(np.exp(-1.*(10**nu_break/nu)**(a1-a2)))) / \
        (1.-np.float128(np.exp(-1.*(10**nu_to/nu)**(a0-a1))))
#    y /= (1+ redshift)
#    y /= nu
    return np.float64(y)

def abs_double_sync_cutoff_law(nu, param, redshift):
    """
    calculate the flux for a double powerlaw, absorption at low frequency and cut-off at 
    high frequency. Contains 3 spectral index and 3 break frequencies.

    :param nu: array of increasing frequency (in log)
    :param param: array of the parameter values
    :param redshift: redshift of the source
    :return: the calculated flux at the provided frequency array with the given parameters

    parameters array description:
    p[0]: norm => normalisation factor in log
    p[1]: nu_to => frequency where tau=1 (optically thick/thin transition) in log
    p[2]: nu_break => break in higher frequency spectral index in log
    p[3]: nu_co => cut-off frequency in log
    p[4]: alpha0 => absorption index
    p[5]: alpha1 => lower frequency index
    p[6]: alpha2 => higher frequency index
    """

    return y

def Polletta2000_SSRQ(nu, param, redshift):
    """
    Calculate the flux for a double power law and cutoff at high frequency.
    Inspired from Polletta et al. 200 formula. Contains 3 spectral slopes and 2 break frequencies

    :param nu: array of increasing frequency (in log)
    :param param: array of the parameter values
    :param redshift: redshift of the source
    :return: the calculated flux at the provided frequency array with the given parameters

    parameters array description:
    p[0]: norm => normalisation factor in log
    p[1]: nu_t => frequency where tau=1 (optically thick/thin transition) in log
    p[2]: nu_cutoff => high frequency cutoff (plasma energy) in log
    p[3]: alpha1 => optically thick index
    p[4]: alpha2 => optically thin index
    """

    n, nut, nuco, a1, a2 = param
    y = np.float128(10**n*(nu*(1.+redshift)/10**nut)**a1) * \
        (1.-np.float128(np.exp(-1.*(10**nut/(nu*(1.+redshift)))**(a1-a2)))) * \
        np.float128(np.exp(-1.0*nu*(1.+redshift)/10**nuco))
    # to check
    y /= (1+redshift)
    y /= nu
    return np.float64(y)


def AGN_law(nu, param, redshift):
    # check redshift
    freq_cut = 3e14 / 25.  # 25 in um
    norm, alpha = param
    nu_z = nu * (1. + redshift)
    y = (freq_cut/nu_z)**(-alpha) * np.exp(-freq_cut/nu_z)*(10**norm)
    y /= (1+ redshift)
    return y


def triple_sync_law(nu,param,redshift):
    """
    Calculate the flux for a triple power law with two break frequencies. Inspired from the formula
    at https://math.stackexchange.com/questions/2427089/how-do-i-smoothly-merge-two-power-laws
    
    :param nu: array of increasing frequency
    :param param: array of the parameter values
    :param redshift: redshift of the source
    :return: the calculated flux at the provided frequency for the given parameter

    parameters array description:
    p[0]: norm => normalisation
    p[1]: nu_b1 => low frequency break
    p[2]: nu_b2 => high frequency break
    p[3]: a1 => lower frequency index
    p[4]: a2 => mid frequency index
    p[5]: a3 => high frequency index
    """
    norm,nu_b1,nu_b2,a1,a2,a3=param
    s=1.
    y = 10**norm * (nu*(1.+redshift))**a1 * \
        (1.+(nu*(1.+redshift)/10**nu_b1)**(np.abs(a1-a2)*s))**(np.sign(a2-a1)/s) * \
        (1.+(nu*(1.+redshift)/10**nu_b2)**(np.abs(a2-a3)*s))**(np.sign(a3-a2)/s)
    y /= (1.+redshift)
    return y


def triple_sync_law2(nu,param,redshift):
    """
    Calculate the flux for a triple power law with two break frequencies. Inspired from the formula
    at https://math.stackexchange.com/questions/2427089/how-do-i-smoothly-merge-two-power-laws
    
    :param nu: array of increasing frequency
    :param param: array of the parameter values
    :param redshift: redshift of the source
    :return: the calculated flux at the provided frequency for the given parameter

    parameters array description:
    p[0]: norm => normalisation
    p[1]: nu_b1 => low frequency break
    p[2]: nu_b2 => high frequency break
    p[3]: a1 => lower frequency index
    p[4]: a2 => mid frequency index
    p[5]: a3 => high frequency index
    """
    norm,nu_b1,nu_b2,a1,a2,a3=param
    s=1.
    # trick to make it more stable during fitting; change of referential point
    ref_point=np.abs(nu_b1-nu_b2)+min(nu_b1,nu_b2)
    nu = nu / 10**ref_point
    y = 10**norm*(nu*(1.+redshift))**a1 * \
        (1.+(nu*(1.+redshift)/10**(nu_b1-ref_point))**(np.abs(a1-a2)*s))**(np.sign(a2-a1)/s) * \
        (1.+(nu*(1.+redshift)/10**(nu_b2-ref_point))**(np.abs(a2-a3)*s))**(np.sign(a3-a2)/s)
    y /= (1.+redshift)
    return y


def modified_BB_draine2007(nu,param,redshift):
    """
    Formula of modified black body to link dust mass
    and beta (and kappa), without an extra normalisation factor
    """
    import astropy.constants as cst
    from astropy.cosmology import Planck15
    from astropy import analytic_functions
    import astropy.units as u
    from astropy.modeling import models
    mdust, tdust = param
    bb=models.BlackBody1D(tdust*u.K,bolometric_flux=cst.L_sun.to(u.erg/u.s)/(u.cm*u.cm))
    nu_rest = nu * (1+redshift) * u.Hz # shift the frequency range to restframe
    nu_0 = (cst.c/(250e-6*u.m)).to(u.Hz) 
    kappa_n0= 4.0*u.cm*u.cm/u.g
    beta=2.08
    kappa_abs = kappa_n0.to(u.m*u.m/u.kg) * (nu_rest/nu_0)**beta
    DL = Planck15.luminosity_distance(redshift)

    y_z = bb(nu_rest) * (1+redshift) * kappa_abs.to(u.cm*u.cm/u.kg)
    y_z *= ((10**mdust)*u.M_sun).to(u.kg)*u.kg
    y_z /= (DL.to(u.cm) * DL.to(u.cm) * cst.M_sun)
    return  y_z.value

def modified_BB_draine2007_z(nu,param):
    """
    Formula of modified black body to link dust mass
    and beta (and kappa), without an extra normalisation factor
    """
    import astropy.constants as cst
    from astropy.cosmology import Planck15
    from astropy import analytic_functions
    import astropy.units as u
    from astropy.modeling import models
    mdust, tdust, redshift = param
    bb=models.BlackBody1D(tdust*u.K,bolometric_flux=cst.L_sun.to(u.erg/u.s)/(u.cm*u.cm))
    nu_rest = nu * (1+redshift) * u.Hz # shift the frequency range to restframe
    nu_0 = (cst.c/(250e-6*u.m)).to(u.Hz) 
    kappa_n0= 4.0*u.cm*u.cm/u.g
    beta=2.08
    kappa_abs = kappa_n0.to(u.m*u.m/u.kg) * (nu_rest/nu_0)**beta
    DL = Planck15.luminosity_distance(redshift)

    y_z = bb(nu_rest) * (1+redshift) * kappa_abs.to(u.cm*u.cm/u.kg)
    y_z *= ((10**mdust)*u.M_sun).to(u.kg)*u.kg
    y_z /= (DL.to(u.cm) * DL.to(u.cm) * cst.M_sun)
    return  y_z.value

def modified_BB_draine2007beta_z(nu,param):
    """
    Formula of modified black body to link dust mass
    and beta (and kappa), without an extra normalisation factor
    """
    import astropy.constants as cst
    from astropy.cosmology import Planck15
    from astropy import analytic_functions
    import astropy.units as u
    from astropy.modeling import models
    mdust, tdust, beta, redshift = param
    bb=models.BlackBody1D(tdust*u.K,bolometric_flux=cst.L_sun.to(u.erg/u.s)/(u.cm*u.cm))
    nu_rest = nu * (1+redshift) * u.Hz # shift the frequency range to restframe
    nu_0 = (cst.c/(250e-6*u.m)).to(u.Hz) 
    kappa_n0= 4.0*u.cm*u.cm/u.g
    kappa_abs = kappa_n0.to(u.m*u.m/u.kg) * (nu_rest/nu_0)**beta
    DL = Planck15.luminosity_distance(redshift)

    y_z = bb(nu_rest) * (1+redshift) * kappa_abs.to(u.cm*u.cm/u.kg)
    y_z *= ((10**mdust)*u.M_sun).to(u.kg)*u.kg
    y_z /= (DL.to(u.cm) * DL.to(u.cm) * cst.M_sun)
    return  y_z.value
