import readsubfHDF5
import h5py
import numpy as np

from units import AREPO_units
AU = AREPO_units()

import matplotlib
matplotlib.use('PDF')
from matplotlib import rc
#rc('text', usetex=True)
rc('font', family='serif')
import matplotlib.pyplot as plt
import matplotlib.ticker
import brewer2mpl

base = "/n/ghernquist/Illustris/Runs/L75n1820FP/"

# Playing around with group catalog
# note the following code only works when logged in on dev
# cat=readsubfHDF5.subfind_catalog("/scratch/sims.illustris/L75n1820FP",135)
cat = readsubfHDF5.subfind_catalog(base,120,keysel=['SubhaloGrNr','GroupFirstSub'])
#135 - z=0
#120 - z=0.19728

# Playing around with galaxy properties catalog
galf = h5py.File(base+"/postprocessing/galprop/galprop_120.hdf5",'r')
galf.keys()
gal_sm = AU.PhysicalMass(np.array(galf['stellar_totmass']))
gal_sfr = np.array(galf['gas_sfr_inrad'])
n_gal = np.size(gal_sm)
sm_mass_select = np.logical_and(gal_sm > 10.**9.5,gal_sm < 10.**11.5)
print np.sum(sm_mass_select) #5274
gal_ids = np.arange(n_gal)[sm_mass_select]

# filter out only primary halos:
# method 1 - must be most massive stellar-wise in its group [not implemented yet]
# method 2 - must inhabit most massive subhalo in the group
primary_gal_ids = np.ones(0)
sat_gal_ids = np.ones(0)
sat_gal_grnr = np.ones(0)
for gal_id in gal_ids:
	grnr = cat.SubhaloGrNr[gal_id]
	if cat.GroupFirstSub[grnr] ==  gal_id: 
		primary_gal_ids = np.append(primary_gal_ids,gal_id)
	else:
		sat_gal_ids = np.append(sat_gal_ids,gal_id)
		sat_gal_grnr = np.append(sat_gal_grnr,grnr)


grp_cent_sfr = np.zeros(0)
grp_sat_sfr = np.zeros(0)
grp_n_sat = np.zeros(0)
grp_mostSFR_satsm = np.zeros(0)
for grnr in cat.SubhaloGrNr[np.int32(primary_gal_ids)]:
	# if grnr < 200:
	this_grp_cent_sfr = gal_sfr[cat.GroupFirstSub[grnr]]
	sat_ids = sat_gal_ids[sat_gal_grnr == grnr]
	# print "sat_ids ",sat_ids
	sat_sfr = gal_sfr[np.int32(sat_ids)]
	sat_sm = gal_sm[np.int32(sat_ids)]
	this_grp_sat_sfr = np.sum(gal_sfr[np.int32(sat_ids)])
	print "Group: {} - Gal SFR {} - Sat SFR {} - Ratio {} ".format(grnr,this_grp_cent_sfr,this_grp_sat_sfr,np.log10(this_grp_sat_sfr/this_grp_cent_sfr))
	grp_cent_sfr = np.append(grp_cent_sfr,this_grp_cent_sfr)
	grp_sat_sfr = np.append(grp_sat_sfr,this_grp_sat_sfr)
	grp_n_sat = np.append(grp_n_sat,np.size(sat_ids))
	if np.size(sat_ids) > 0:
		grp_mostSFR_satsm = np.append(grp_mostSFR_satsm,sat_sm[np.argmax(sat_sfr)])
	else:
		grp_mostSFR_satsm = np.append(grp_mostSFR_satsm,0.)

x = np.log10(gal_sm[np.int32(primary_gal_ids)])
y = np.log10(grp_sat_sfr/grp_cent_sfr)
y[grp_sat_sfr == 0.] = -4.
c = np.zeros_like(x)
ms = np.log10(grp_mostSFR_satsm)
c[np.logical_and(ms > 11.,ms < 12.)] = 1.
c[np.logical_and(ms > 10.,ms < 11.)] = 2.
c[np.logical_and(ms > 9.,ms < 10.)] = 3.
# c[np.logical_and(ms > 8.,ms < 9.)] = 4.
# c[np.logical_and(ms > 7.,ms < 8.)] = 5.
# bmap = brewer2mpl.get_map('Spectral','Diverging',9)
# cmap = bmap.get_mpl_colormap(N=5, gamma=2.0)
# c = cmap(c/4.)
c *= 2.

plt.scatter(x,y,c=c,marker='+',alpha=0.6)
plt.xlabel('Log10[Central Gal Mass]')
plt.ylabel('Log10[Total Satellite SFR/Central SFR]')
plt.savefig('satellite_SFR.pdf')

plt.close('all')
x = np.log10(gal_sm[np.int32(primary_gal_ids)])
y = grp_n_sat
plt.scatter(x,y)
plt.xlabel('Log10[Central Gal Mass]')
plt.ylabel('# of Satellites')
plt.savefig('satellite_num.pdf')
