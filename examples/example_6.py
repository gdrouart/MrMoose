"""
Example to generate a .fit, .mod and .dat file to feed in MrMoose for 
demonstration. The model consists of two single power-laws, two black bodies and
a AGN modified black body, with 18 data points of which 3 are an upper limits, 
and a mixture of unresolved and blended/spatially identified components, all at z=2. 
This example is inspired from a real case from Falkendal et al. 2018 paper. 
"""

import models as md
import numpy as np
import mm_utilities as mm
import read_files as rd

norm_sync1 = 1.0  
alpha_sync1 = -1.
norm_sync2 = 1.0
alpha_sync2 = -1.2
norm_bb1 = 1.0
temp1 = 40
norm_bb2 = 1.0
temp2 = 60
norm_agn = -9.
alpha_agn = -2.

nu = 10**np.linspace(6, 18, 10000)  # frequency range
redshift = [2.0, 2.0, 2.0, 2.0, 2.0]

# generate with the provided model
fnu = md.sync_law(nu, [norm_sync1, alpha_sync1], redshift[0]) + \
      md.sync_law(nu, [norm_sync2, alpha_sync2], redshift[1]) + \
      md.BB_law(nu, [norm_bb1, temp1], redshift[2]) + \
      md.BB_law(nu, [norm_bb2, temp2], redshift[3]) + \
      md.AGN_law(nu, [norm_agn, alpha_agn], redshift[4])

# list of the filters
filter_name = np.array(['VLA_L', 'VLA_C', 'VLA_C', 'VLA_X', 'VLA_X',
                        'ATCA_47', 'ALMA_3', 'ALMA_6', 'ALMA_6_nr1',
                        'laboca_870', 'spire_500', 'spire_350', 'spire_250',
                        'pacs_160', 'pacs_70', 'mips_24', 'irs_16', 'irac_4'])
data_nature = np.array(['d', 'd', 'd', 'd', 'd',
                        'd', 'd', 'd', 'd',
                        'u', 'u', 'u', 'd',
                        'd', 'd', 'd', 'd', 'd'])  # "d" for detections, "u" for upper limit
arrangement = np.array(['4', '6', '5', '6', '5',
                        '4', '7', '3', '2',
                        '1', '1', '1', '1',
                        '1', '1', '1', '1', '1'])  # do not forget the "," for the last element!
comp_number = np.array(['3,4',   '4',         '3',     '4',      '3',
                        '3,4',   '0,1,2,3,4', '2,4',   '1,3',
                        '0,1,2', '0,1,2',     '0,1,2', '0,1,2',
                        '0,1,2', '0,1,2',     '0,1,2', '0,1,2', '0,1,2'])
sn_mod = np.array([5., 5., 5., 5., 5.,
                   5., 5., 5., 5.,
                   5., 5., 5., 5.,
                   5., 5., 5., 5., 5.])  # SN detection to estimate noise level for each point
notes = np.array(["'total sync'", "'E sync'", "'W sync'", "'E sync'", "'W sync'",
                  "'total sync'", "'all sources'", "'E sync/host'", "'W sync/comp'",
                  "'Dust/AGN'", "'Dust/AGN'", "'Dust/AGN'", "'Dust/AGN'",
                  "'Dust/AGN'", "'Dust/AGN'", "'Dust/AGN'", "'Dust/AGN'", "'Dust/AGN'"])  # notes on observations
RA_list = ['12h00m00s', '12h00m00.1s', '11h59m59.95s', '12h00m00.1s', '11h59m59.95s',
           '12h00m00s', '12h00m00s', '12h00m00.1s', '11h59m59.95s',
           '12h00m00s', '12h00m00s', '12h00m00s', '12h00m00s',
           '12h00m00s', '12h00m00s', '12h00m00s', '12h00m00s', '12h00m00s']
Dec_list = ['-40d00m00s', '-39d59m59s', '-40d00m01s', '-39d59m59s', '-40d00m01s',
            '-40d00m00s', '-40d00m00s', '-39d59m59s', '-40d00m00.5s',
            '-40d00m00s', '-40d00m00s', '-40d00m00s', '-40d00m00s',
            '-40d00m00s', '-40d00m00s', '-40d00m00s', '-40d00m00s', '-40d00m00s']
res_list = [20., 1.0, 1.0, 0.5, 0.5,
            10., 3.0, 0.3, 0.3,
            15., 35., 25., 17.,
            5., 4., 7.0, 5.0, 5.0]

# create the array to feed in the data file
fnu_mod = np.zeros(filter_name.size)
fnu_err = np.zeros(filter_name.size)
lambda0 = np.zeros(filter_name.size)

# run through the filters to create the simulated data
for i_filter, name_filter in enumerate(filter_name):
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
with open('data/fake_source_ex6.dat', 'w') as fake:
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
with open('fake_source_ex6.fit', 'w') as fake:
    fake.write('source_file: data/fake_source_ex6.dat \n')
    fake.write('model_file: models/fake_source_ex6.mod \n')
    fake.write('all_same_redshift: True \n')
    fake.write('redshift: '+str(redshift)+'\n')
    fake.write('nwalkers: 40 \n')
    fake.write('nsteps: 20 \n')
    fake.write('nsteps_cut: 19 \n')
    fake.write('percentiles: [10., 25., 50., 75., 90.] \n')
    fake.write('skip_imaging: False \n')
    fake.write('skip_fit: False \n')
    fake.write('skip_MCChains: False \n')
    fake.write('skip_triangle: False \n')
    fake.write('skip_SED: False \n')
    fake.write("unit_obs: 'Hz' \n")
    fake.write("unit_flux: 'Jy' \n")

# create the model file
with open('models/fake_source_ex6.mod', 'w') as fake:
    fake.write('AGN_law  2 \n')
    fake.write('$N_{AGN}$   -38  -28 \n')
    fake.write('$\\alpha_{AGN}$ -4.0  0.0 \n')
    fake.write('BB_law  2 \n')
    fake.write('$N_{BB1}$   -28  -18 \n')
    fake.write('$Temp_1$  10  60 \n')
    fake.write('BB_law  2 \n')
    fake.write('$N_{BB2}$   -28  -18 \n')
    fake.write('$Temp_2$  10  60 \n')
    fake.write('sync_law  2 \n')
    fake.write('$N_{s1}$   -28  -18 \n')
    fake.write('$\\alpha_{s1}$ -2.0  0.0 \n')
    fake.write('sync_law  2 \n')
    fake.write('$N_{s2}$   -28  -18 \n')
    fake.write('$\\alpha_{s2}$ -2.0  0.0 \n')
