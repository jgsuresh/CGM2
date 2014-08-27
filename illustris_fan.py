

class illustris_fan
	def __init__(self):
		self.snapbase = "/n/ghernquist/Illustris/Runs/Illustris-1/output/"
		self.CGMsnap_base = '/n/home04/jsuresh/scratch1/AREPOfest/data/CGM_snaps/'
		self.grid_base = '/n/home04/jsuresh/scratch1/AREPOfest/data/grids/'
		self.fig_base = '/n/home04/jsuresh/scratch1/AREPOfest/data/figs/'

		self.snapnum = 120

		#read in catalog
        self.cat = readsubfHDF5.subfind_catalog(self.snapbase, self.snapnum, long_ids=long_ids, double_output=double_output, keysel=["GroupFirstSub","SubhaloGrNr"])
		self.load_gal_props()
		self.gal_mass_vs_sSFR()


	def gal_mass_vs_sSFR(self,include_COS=True,savename=None):
		plt.figure()

		ax = plt.scatter(self.gal_sm,self.gal_ssfr,marker='.')

		if include_COS:
			werk_dat = np.load_txt(werk.dat)
			COS_log_sm = werk_dat[:,0]
			COS_sfr = werk_dat[:,1]
			upper_lim = np.array(werk_dat[:,2],dtype='bool')
			not_upper_lim = np.logical_not(upper_lim)

			COS_sm = 10.**COS_log_sm
			COS_ssfr = COS_sfr/COS_sm
			plt.scatter(COS_sm[not_upper_lim],COS_ssfr[not_upper_lim],marker='o',color='purple')
			plt.scatter(COS_sm[upper_lim],COS_ssfr[upper_lim],marker='v',color='purple')

		ax.set_xscale('log')
		ax.set_yscale('log')
		if savename != None: 
            plt.savefig(self.fig_base+savename+".pdf", bbox_inches='tight')
        else:
            plt.savefig(self.fig_base+"sm_vs_ssfr.pdf", bbox_inches='tight')



	def load_gal_props(self):
		galf = h5py.File(base+"/postprocessing/galprop/galprop_120.hdf5",'r')
		gal_sm = np.array(galf['stellar_totmass'])
		gal_sfr = np.array(galf['gas_totsfr'])

		sm_mass_select = np.logical_and(gal_sm > 1.,gal_sm < 40)
		gal_ids = np.arange(n_gal)[sm_mass_select]

		# filter out only primary halos:
		# method 1 - must be most massive stellar-wise in its group [not implemented yet]
		# method 2 - must inhabit most massive subhalo in the group
		primary_gal_ids = np.ones(0)
		for gal_id in gal_ids:
			grnr = self.cat.SubhaloGrNr[gal_id]
			if self.cat.GroupFirstSub[grnr] ==  gal_id: 
				primary_gal_ids = np.append(primary_gal_ids,gal_id)
		print "np.size(primary_gal_ids) ",np.size(primary_gal_ids)

		# Sanity check: ensure there are no duplicates in the group number associated with the primary_gal_ids
		primary_gal_ids = np.int32(primary_gal_ids)
		grnr = cat.SubhaloGrNr[primary_gal_ids]
		if np.size(np.unique(grnr)) < np.size(grnr): print "had non-uniques!"

		self.gal_sm = AU.PhysicalMass(gal_sm[primary_gal_ids])
		self.gal_sfr = gal_sfr[primary_gal_ids]
		self.gal_ssfr = self.gal_sfr/self.gal_sm


if __name__ == '__main__':
    illustris_fan()
