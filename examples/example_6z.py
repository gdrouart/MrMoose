"""
Example to generate a .fit, .mod and .dat file to feed in MrMoose for 
demonstration. The model consists of one single power-laws and two black 
bodies, with 15 data points. All is a mixture of unresolved and 
blended/spatially identified components, with the black bodies being at 
different redshifts (z=2 and z=4).
"""

from models import *
import numpy as np
import mm_utilities as mm
import read_files as rd

# first group of component at same redshift
redshift1 = 2.0
func1a = 'sync_law'  # need to make sure function is existing in model.py
norm_sync1 = 7.0  # parameters - normalisation
alpha_sync1 = -2.  # parameters - spectral index
func1b = 'BB_law'
norm_bb1 = 1.0  # parameter - normalisation
temp1 = 40  # parameter - temperature [K]

# second group of component at other redshift
redshift2 = 4.0
func2a = 'BB_law'
norm_bb2 = 0.1  # parameter - normalisation
temp2 = 20  # parameter - temperature [K]

# making all in form to build the fake system
# array of the function name
comp_function = np.array([func1a, func1b, func2a])
# array of the redshift of the component
comp_redshift = np.array([redshift1, redshift1, redshift2])
# array of parameter values, organised as sub-arrays respecting function calls
comp_param = np.array([[norm_sync1, alpha_sync1], [norm_bb1, temp1], [norm_bb2, temp2]])

nu = 10**np.linspace(6, 18, 10000)  # frequency range

# list of the filters, arrangements and components
filter_name = np.array(['VLA_L', 'VLA_C', 'VLA_C', 'VLA_X', 'VLA_X',
                        'ATCA_47', 'ALMA_3', 'ALMA_6', 'ALMA_6_nr1',
                        'laboca_870', 'spire_500', 'spire_350', 'spire_250',
                        'pacs_160', 'pacs_70'])
data_nature = np.array(['d', 'd', 'd', 'd', 'd',
                        'd', 'd', 'd', 'd',
                        'd', 'd', 'd', 'd',
                        'd', 'd'])  # "d" for detections, "u" for upper limit
arrangement = np.array(['1', '1', '1', '1', '1',
                        '1', '2', '3', '4',
                        '5', '5', '5', '5',
                        '5', '5'])  # do not forget the "," for the last element!
comp_number = np.array(['0', '0', '0', '0', '0',
                        '0', '0,1,2', '0,1', '2',
                        '1,2', '1,2', '1,2', '1,2',
                        '1,2', '1,2'])
sn_mod = np.array([5., 5., 5., 5., 5.,
                   5., 5., 5., 5.,
                   5., 5., 5., 5.,
                   5., 5.])  # SN detection to estimate noise level for each point
notes = np.array(["'sync'", "'sync'", "'sync'", "'sync'", "'sync'",
                  "'sync'", "'all'", "'host'", "'comp'",
                  "'host+comp'", "'host+comp'", "'host+comp'", "'host+comp'",
                  "'host+comp'", "'host+comp'"])  # notes on observations
RA_list = ['12h00m00s', '12h00m00.1s', '11h59m59.95s', '12h00m00.1s', '11h59m59.95s',
           '12h00m00s', '12h00m00s', '12h00m00.1s', '11h59m59.95s',
           '12h00m00s', '12h00m00s', '12h00m00s', '12h00m00s',
           '12h00m00s', '12h00m00s']
Dec_list = ['-40d00m00s', '-39d59m59s', '-40d00m01s', '-39d59m59s', '-40d00m01s',
            '-40d00m00s', '-40d00m00s', '-39d59m59s', '-40d00m00.5s',
            '-40d00m00s', '-40d00m00s', '-40d00m00s', '-40d00m00s',
            '-40d00m00s', '-40d00m00s']
res_list = [20., 1.0, 1.0, 0.5, 0.5,
            10., 3.0, 0.3, 0.3,
            15., 35., 25., 17.,
            5., 4.]

# create the array to feed in the data file
fnu_mod = np.zeros(filter_name.size)
fnu_err = np.zeros(filter_name.size)
lambda0 = np.zeros(filter_name.size)

# convert the component numbers into integer list to create the combined SED following the provided arrangements
func_index = [map(int, (elem.replace(',', ''))) for elem in comp_number]

# run through the filters to create the simulated data
for i_filter, name_filter in enumerate(filter_name):
    # calculate the sum of components for this arrangement
    fnu = [globals()[comp_function[j]](nu, comp_param[j], comp_redshift[j]) for j in func_index[i_filter]]
    # trick to get rid off the extra dimension
    fnu = np.sum(fnu, axis=0)

    # read the filter transmission
    nu_filter, trans_filter = rd.read_single_filter('filters/'+name_filter+'.fil')
    # calculate the lambda0
    lambda0[i_filter] = np.average(nu_filter, weights=trans_filter)
    # perform the integration
    tmp = mm.integrate_filter(nu, fnu, nu_filter, trans_filter)
    # add a gaussian noise (depending on the signal to noise defined previously)
    fnu_err[i_filter] = tmp/sn_mod[i_filter]
    fnu_mod[i_filter] = np.random.normal(tmp, fnu_err[i_filter])
    if data_nature[i_filter] == 'u':
        fnu_err[i_filter] = fnu_mod[i_filter]


# create the data file
with open('data/fake_source_ex6z.dat', 'wb') as fake:
    fake.writelines("# filter        RA              Dec        resolution  lambda0  det_type  flux   "
                    "flux_error  arrangement  component   component_number \n")
    for i in range(filter_name.size-1):
        fake.write('{:15} {:15} {:15} {:5.1f} {:10e} {:5} {:10e} {:10e} {:10} {:10} {:10} \n'.format(
            filter_name[i], RA_list[i], Dec_list[i], res_list[i],
            lambda0[i], data_nature[i], fnu_mod[i], fnu_err[i], arrangement[i], notes[i], comp_number[i]))
    fake.write('{:15} {:15} {:15} {:5.1f} {:10e} {:5} {:10e} {:10e} {:10} {:10} {:10}'.format(
        filter_name[i+1], RA_list[i+1], Dec_list[i+1], res_list[i+1],
        lambda0[i+1], data_nature[i+1], fnu_mod[i+1], fnu_err[i+1], arrangement[i+1], notes[i+1], comp_number[i+1]))

# create the fit file
with open('fake_source_ex6z.fit', 'wb') as fake:
    fake.write('source_file: data/fake_source_ex6z.dat \n')
    fake.write('model_file: models/fake_source_ex6z.mod \n')
    fake.write('all_same_redshift: False \n')
    fake.write('redshift: '+"[{:.4f}, {:.4f}, {:.4f}]".format(redshift1, redshift1, redshift2)+'\n')
    fake.write('nwalkers: 20 \n')
    fake.write('nsteps: 80 \n')
    fake.write('nsteps_cut: 78 \n')
    fake.write('percentiles: [10., 25., 50., 75., 90.] \n')
    fake.write('skip_imaging: False \n')
    fake.write('skip_fit: False \n')
    fake.write('skip_MCChains: False \n')
    fake.write('skip_triangle: False \n')
    fake.write('skip_SED: False \n')
    fake.write("unit_obs: 'Hz' \n")
    fake.write("unit_flux: 'Jy' \n")

# create the model file
with open('models/fake_source_ex6z.mod', 'wb') as fake:
    fake.write('sync_law  2 \n')
    fake.write('$N_{s1}$   -22  -12 \n')
    fake.write('$\\alpha_{s1}$ -3.5  -0.5 \n')
    fake.write('BB_law  2 \n')
    fake.write('$N_{BB1}$   -28  -18 \n')
    fake.write('$T_1$  10  60 \n')
    fake.write('BB_law  2 \n')
    fake.write('$N_{BB2}$   -28  -18 \n')
    fake.write('$T_2$  10  40 \n')
