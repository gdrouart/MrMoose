import models as md
import numpy as np
import mm_utilities as mm
import read_files as rd

#def fake_sync_source():
# define the parameters of the model and create
# the normalisation will effectively be  values-23 due to the Jansky transformation
norm_sync = 1.0
alpha_sync = -1.0
norm_bb = 1.0
temp = 40
norm_agn = -9.
alpha_agn = -2.

nu = 10**np.linspace(6, 18, 10000)  # frequency range
redshift = 0

# generate with the provided model
fnu = md.sync_law(nu, [norm_sync, alpha_sync], redshift) + \
      md.BB_law(nu, [norm_bb, temp], redshift) + \
      md.AGN_law(nu, [norm_agn, alpha_agn], redshift)
#print md.sync_law(nu, [norm_sync, alpha], redshift)
#print
#print md.BB_law(nu, [norm_bb, temp], redshift)

# list of the filters
filter_name = np.array(['74MHz(VLA)', '408MHz',
                        'ALMA_3', 'ALMA_6', 'laboca_870', 'spire_500', 'spire_350',
                        'spire_250', 'pacs_160', 'pacs_100', 'pacs_70', 'mips_24',
                        'irac_4', 'irac_3', 'irac_2', 'irac_1'])
data_nature = np.array(['d', 'd',
                        'u', 'd', 'd', 'd', 'd',
                        'd', 'd', 'd', 'd', 'd',
                        'd', 'd', 'd', 'd'])  # "d" for detections, "ul" for upper limit
arrangement = np.array(['1', '1',
                        '1', '1', '1', '1', '1',
                        '1', '1', '1', '1', '1',
                        '1', '1', '1', '1'])  # do not forget the "," for the last element!
comp_number = np.array(['0,1,2', '0,1,2',
                        '0,1,2', '0,1,2', '0,1,2', '0,1,2', '0,1,2',
                        '0,1,2', '0,1,2', '0,1,2', '0,1,2', '0,1,2',
                        '0,1,2', '0,1,2', '0,1,2', '0,1,2'])
sn_mod = [15., 15.,
          15., 15., 15., 15., 15.,
          15., 15., 15., 15., 15.,
          15., 15., 15., 15.]  # SN detection to estimate noise level for each point

RA_list = ['12h00m00s', ]*16
Dec_list = ['-40d00m00s', ]*16
res_list = [50., 35.,
            0.5, 0.5, 15., 35., 25.,
            15., 5., 4., 4., 7.,
            4., 3.5, 3.0, 3.0]

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
with open('data/fake_source_ex4.dat', 'wb') as fake:
    fake.writelines("# filter        RA              Dec        resolution  lambda0  det_type  flux   "
                    "flux_error  arrangement  component   component_number \n")
    for i in range(filter_name.size-1):
        fake.write('{:15} {:15} {:15} {:5.1f} {:10e} {:5} {:10e} {:10e} {:10} {:10} {:10} \n'.format(
            filter_name[i], RA_list[i], Dec_list[i], res_list[i],
            lambda0[i], data_nature[i], fnu_mod[i], fnu_err[i], arrangement[i], "note", comp_number[i]))
    fake.write('{:15} {:15} {:15} {:5.1f} {:10e} {:5} {:10e} {:10e} {:10} {:10} {:10}'.format(
        filter_name[i+1], RA_list[i+1], Dec_list[i+1], res_list[i+1],
        lambda0[i+1], data_nature[i+1], fnu_mod[i+1], fnu_err[i+1], arrangement[i+1], "note", comp_number[i+1]))

# create the fit file
with open('fake_source_ex4.fit', 'wb') as fake:
    fake.write('source_file: data/fake_source_ex4.dat \n')
    fake.write('model_file: models/fake_source_ex4.mod \n')
    fake.write('redshift: '+str(redshift)+'\n')
    fake.write('nwalkers: 100 \n')
    fake.write('nsteps: 400 \n')
    fake.write('nsteps_cut: 380 \n')
    fake.write('percentiles: [10., 25., 50., 75., 90.] \n')
    fake.write('skip_imaging: False \n')
    fake.write('skip_fit: False \n')
    fake.write('skip_MCChains: False \n')
    fake.write('skip_triangle: False \n')
    fake.write('skip_SED: False \n')
    fake.write("unit_obs: 'Hz' \n")
    fake.write("unit_flux: 'Jy' \n")

# create the model file
with open('models/fake_source_ex4.mod', 'wb') as fake:
    fake.write('sync_law  2 \n')
    fake.write('$N_s$   -28  -18 \n')
    fake.write('$\\alpha_s$ -2.0  0.0 \n')
    fake.write('BB_law  2 \n')
    fake.write('$N_{BB}$   -28  -18 \n')
    fake.write('$Temp$  10  60 \n')
    fake.write('AGN_law  2 \n')
    fake.write('$N_{AGN}$   -38  -28 \n')
    fake.write('$\\alpha_{AGN}$ -4.0  0.0 \n')
