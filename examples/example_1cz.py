"""
Example to generate a .fit, .mod and .dat file to feed in MrMoose for 
demonstration. The model consists of a double power-law and and a break frequency
six data points from a source at z=4, but will be fitted as a free parameter
"""
import sys
# adding the path 
sys.path.insert(0, '/Users/guillaume/Desktop/MrMoose/MrMoose/')

import utils.models as md
import numpy as np
import utils.mm_utilities as mm
import utils.read_files as rd

norm = 1.0
nu_break = 9.0
alpha1 = 2.
alpha2 = -1.5

nu = 10**np.linspace(6, 11, 10000)
redshift = 4.
fnu = md.double_sync_law(nu, [norm, nu_break, alpha1, alpha2], redshift)

filter_name = np.array(['74MHz(VLA)', '178MHz', '408MHz', '1.4GHz', '4.85GHz', '8.4GHz'])
sn_mod = [15., 15., 15., 15., 15., 15.]
RA_list = ['12h00m00s', ]*6
Dec_list = ['-40d00m00s', ]*6
res_list = [12., ]*6

fnu_mod = np.zeros(filter_name.size)
fnu_err = np.zeros(filter_name.size)
lambda0 = np.zeros(filter_name.size)

# run through the filters
for i_filter, name_filter in enumerate(filter_name):
    # read the filter transmission
    nu_filter, trans_filter = rd.read_single_filter('../filters/'+name_filter+'.fil')
    # calculate the lambda0
    lambda0[i_filter] = np.average(nu_filter, weights=trans_filter)
    # perform the integration
    tmp = mm.integrate_filter(nu, fnu, nu_filter, trans_filter)
    # add a gaussian noise (depending on the signal to noise defined previously)
    fnu_err[i_filter] = tmp/sn_mod[i_filter]
    fnu_mod[i_filter] = np.random.normal(tmp, fnu_err[i_filter])

# create the data file
with open('../data/fake_source_ex1cz.dat', 'w') as fake:
    fake.write("# filter        RA              Dec        resolution  lambda0  det_type  flux   flux_error  arrangement  component   component_number \n")
    for i_filter in range(filter_name.size-1):
        fake.write('{:15} {:15} {:15} {:5.1f} {:10e} {:5} {:10e} {:10e} {:10} {:10} {:10} \n'.format(
            filter_name[i_filter], RA_list[i_filter], Dec_list[i_filter], res_list[i_filter],
            lambda0[i_filter], "d", fnu_mod[i_filter], fnu_err[i_filter], "1", "note", "0"))
    fake.write('{:15} {:15} {:15} {:5.1f} {:10e} {:5} {:10e} {:10e} {:10} {:10} {:10}'.format(
        filter_name[i_filter+1], RA_list[i_filter+1], Dec_list[i_filter+1], res_list[i_filter+1],
        lambda0[i_filter+1], "d", fnu_mod[i_filter+1], fnu_err[i_filter+1], "1", "note", "0,"))

# create the fit file
redshift_fit = -1.
with open('../fake_source_ex1cz.fit', 'w') as fake:
    fake.write('source_file: data/fake_source_ex1cz.dat \n')
    fake.write('model_file: models/fake_source_ex1cz.mod \n')
    fake.write('all_same_redshift: True \n')
    fake.write('redshift: ['+str(redshift_fit)+,'] \n')
    fake.write('nwalkers: 20 \n')
    fake.write('nsteps: 20 \n')
    fake.write('nsteps_cut: 18 \n')
    fake.write('percentiles: [10., 25., 50., 75., 90.] \n')
    fake.write('skip_imaging: True \n')
    fake.write('skip_fit: False \n')
    fake.write('skip_MCChains: False \n')
    fake.write('skip_triangle: False \n')
    fake.write('skip_SED: False \n')
    fake.write("unit_obs: 'Hz' \n")
    fake.write("unit_flux: 'Jy' \n")

# create the model file
with open('../models/fake_source_ex1cz.mod', 'w') as fake:
    fake.write('double_sync_law_z  5 \n')
    fake.write('$N$       -25  -15 \n')
    fake.write('$\\nu_{break}$  8.0  10.0 \n')
    fake.write('$\\alpha_1$   0.5  4.0 \n')
    fake.write('$\\alpha_2$  -3.0  0.0 \n')
    fake.write('$z$             0  7.0 \n')
