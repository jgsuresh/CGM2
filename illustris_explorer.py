import readsubfHDF5
import h5py
import numpy as np

from units import AREPO_units
AU = AREPO_units()


base = "/n/ghernquist/Illustris/Runs/L75n1820FP/"

# Playing around with group catalog
# note the following code only works when logged in on dev
# cat=readsubfHDF5.subfind_catalog("/scratch/sims.illustris/L75n1820FP",135)
cat = readsubfHDF5.subfind_catalog(base,120)
#135 - z=0
#120 - z=0.19728

# Playing around with galaxy properties catalog
galf = h5py.File(base+"/postprocessing/galprop/galprop_120.hdf5",'r')
galf.keys()
gal_sm = AU.PhysicalMass(np.array(galf['stellar_totmass']))
n_gal = np.size(gal_sm)
sm_mass_select = np.logical_and(gal_sm > 10.**9.5,gal_sm < 10.**11.5)
print np.sum(sm_mass_select) #5274
gal_ids = np.arange(n_gal)[sm_mass_select]

# filter out only primary halos:
# method 1 - must be most massive stellar-wise in its group [not implemented yet]
# method 2 - must inhabit most massive subhalo in the group
primary_gal_ids = np.ones(0)
for gal_id in gal_ids:
	grnr = cat.SubhaloGrNr[gal_id]
	if cat.GroupFirstSub[grnr] ==  gal_id: 
		primary_gal_ids = np.append(primary_gal_ids,gal_id)
print "np.size(primary_gal_ids) ",np.size(primary_gal_ids) #3527 for z=0, 3287 for z=0.2

# Sanity check: ensure there are no duplicates in the group number associated with the primary_gal_ids
primary_gal_ids = np.int32(primary_gal_ids)
grnr = cat.SubhaloGrNr[primary_gal_ids]
if np.size(np.unique(grnr)) < np.size(grnr): print "had non-uniques!"

# Check out corresponding group properties:
print "grnr ",grnr
print "np.max(cat.Group_M_Crit200[grnr]) ",np.max(cat.Group_M_Crit200[grnr]) # 2261.3632812
print "np.min(cat.Group_M_Crit200[grnr]) ",np.min(cat.Group_M_Crit200[grnr]) #7.72121906281
print "np.max(cat.GroupLenType[grnr,4]) ",np.max(cat.GroupLenType[grnr,4]) # 1758605
print "np.min(cat.GroupLenType[grnr,4]) ",np.min(cat.GroupLenType[grnr,4]) # 1758605

# how many pathological cases do we have, where galaxy stellar mass is > 25% of mass of halo?
print "galaxy stellar mass for smallest halo: ", gal_sm[primary_gal_ids[np.argmin(cat.Group_M_Crit200[grnr])]]
print np.sum(gal_sm[primary_gal_ids] > 0.25*cat.Group_M_Crit200[grnr]) # only 1! get rid of this weirdo!

# Find out which subfiles these groups can be found on 


#now find host group for each of the galaxy ids:

# Playing around with individual snapshots
#subf=h5py.File(base+"snapdir_135/snap_135.0.hdf5",'r')
#list(subf['Header'].attrs)
#subf['PartType4'].

