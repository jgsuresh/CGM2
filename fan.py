import os
import time
import numpy as np
import h5py

import matplotlib
matplotlib.use('PDF')
from matplotlib import rc
#rc('text', usetex=True)
rc('font', family='serif')
import matplotlib.pyplot as plt
import matplotlib.ticker
import brewer2mpl
from itertools import cycle

from units import AREPO_units
AU = AREPO_units()
import readsubfHDF5

# import SphMap as sm 

import glob


class full_analysis:
    def __init__(self):
        self.CGMsnap_base = '/n/home04/jsuresh/CGM_new/data/CGM_snaps/'
        self.grid_base = '/n/home04/jsuresh/CGM_new/data/grids/'
        self.fig_base = '/n/home04/jsuresh/CGM_new/data/figs/'
        self.npz_base = '/n/home04/jsuresh/CGM_new/data/npz/'

        # Select which runs to be analyzed concurrently:
        self.run_list = []
        self.color_list = [] 
        self.label_list = [] 
        self.linestyle_list = []
        self.snapnum_list = [] 
        self.snapdir_list = []


        # Feedback paper runs:
        redshifts = ['2']
        self.dat_prep(redshifts,g0_BH=1)
        self.dat_prep(redshifts,c2_256=1)
        self.dat_prep(redshifts,c0_sw_256=1)
        self.dat_prep(redshifts,g50_fixv_nothermal=1)
        self.dat_prep(redshifts,g50_BH=1)
        self.dat_prep(redshifts,g50_fixv=1)
        self.dat_prep(redshifts,c0_nometalwinds=1)
        self.dat_prep(redshifts,c0_fullmetalwinds=1)


        # Check runs:
        # self.dat_prep(redshifts,c0_128=1)
        # self.dat_prep(redshifts,c0_512=1)
        # self.dat_prep(redshifts,g0_BH=1)
        # self.dat_prep(redshifts,c0_check=1)
        # self.dat_prep(redshifts,c0_nometalwinds=1)
        # self.dat_prep(redshifts,c0_fullmetalwinds=1)
        # self.dat_prep(redshifts,c0_dumpfactor95=1)



        # # Cosmo0_V6 (Fiducial):
        # redshifts = ['2']
        # self.dat_prep(redshifts,c0_256=1)
        # self.dat_prep(redshifts,c0_512=1)

        # # Cosmo0_V6_fastWinds, Cosmo0_V6_strongWinds (256^3 particles):
        # redshifts = ['2']
        # self.dat_prep(redshifts,c0_fw_256=1)
        # self.dat_prep(redshifts,c0_sw_256=1)

        # # Cosmo2_V6 (256^3 particles):
        # redshifts = ['2']
        # self.dat_prep(redshifts,c2_256=1)
        
        # # Cosmo3_V6 (No feedback - 512^3 particles):
        # redshifts = ['4','3','2']
        # self.dat_prep(redshifts,c3_512=1)

        # # Cosmo4_V6 (Simeon's heated winds model - 512^3 particles):
        # redshifts = ['2.5']
        # self.dat_prep(redshifts,c4_512=1)
        # self.dat_prep(redshifts,c4_check=1)


        # Gamma boxes (256^3 particles):
        # redshifts = ['2']
        # self.dat_prep(redshifts,g0_BH=1) #c2_256
        # self.dat_prep(redshifts,g10_BH=1)
        # self.dat_prep(redshifts,g20_BH=1)
        # self.dat_prep(redshifts,g30_BH=1)
        # self.dat_prep(redshifts,g40_BH=1)
        # self.dat_prep(redshifts,g50_BH=1)
        # self.dat_prep(redshifts,g10_noBH=1)
        # self.dat_prep(redshifts,g20_noBH=1)
        # self.dat_prep(redshifts,g30_noBH=1)
        # self.dat_prep(redshifts,g40_noBH=1)
        # self.dat_prep(redshifts,g50_noBH=1)

        # redshifts = ['2']
        # self.dat_prep(redshifts,g25_BH=1)
        # self.dat_prep(redshifts,g75_BH=1)
        # self.dat_prep(redshifts,g95_BH=1)
        # self.dat_prep(redshifts,g25_noBH=1)
        # self.dat_prep(redshifts,g75_noBH=1)
        # self.dat_prep(redshifts,g95_noBH=1)

        # redshifts = ['2']
        # self.dat_prep(redshifts,g10_nothermal=1)
        # self.dat_prep(redshifts,g20_nothermal=1)
        # self.dat_prep(redshifts,g30_nothermal=1)
        # self.dat_prep(redshifts,g40_nothermal=1)
        # self.dat_prep(redshifts,g50_nothermal=1)

        # redshifts = ['2']
        # self.dat_prep(redshifts,g25_fixv=1)
        # self.dat_prep(redshifts,g50_fixv=1)

        # redshifts = ['2']
        # self.dat_prep(redshifts,g25_fixv_fixeta=1)
        # self.dat_prep(redshifts,g50_fixv_fixeta=1)

        # redshifts = ['2']
        # self.dat_prep(redshifts,g25_fixv_nothermal=1)
        # self.dat_prep(redshifts,g50_fixv_nothermal=1)
        
        # redshifts = ['2']
        # self.dat_prep(redshifts,c0_nometalwinds=1)
        # self.dat_prep(redshifts,c0_fullmetalwinds=1)


        print self.run_list
        print self.snapnum_list
        print self.color_list
        print self.linestyle_list


        #########################################
        # Full-snapshot data analysis functions #
        #########################################

        ##############################
        # 3D data analysis functions #
        ##############################
        # self.phase_budget('m',savename='m_budget_lowmetal',minimal=False)
        # self.phase_budget('z',savename='z_budget_lowmetal',minimal=False)
        # self.CGM_and_gal_metallicity(savename="massmet_test")
        # self.radial_profile('T')
        # self.radial_profile('z')

        ##############################
        # 2D data analysis functions #
        ##############################
        print "1"
        self.grid_sightlines("H1",200.,coldens_min=17.2,coldens_max=20.3,minmass=10**11.8,maxmass=10**12.2,coverfrac_within_R=True,rudie_LLS=True,savename="LLS_newfan")
        print "2"
        self.grid_sightlines("H1",200.,coldens_min=20.3,minmass=10**11.8,maxmass=10**12.2,coverfrac_within_R=True,rudie_DLA=True,savename="DLA_newfan")
        print "3"
        self.grid_sightlines("H1",200.,minmass=10**11.8,maxmass=10**12.2,savename="coldens_newfan")

    ############################################################################################################
    ############################################################################################################

    #########################################
    # Full-snapshot data analysis functions #
    #########################################

    ##############################
    # 3D data analysis functions #
    ##############################

    ##############################
    # 2D data analysis functions #
    ##############################

    def grid_sightlines(self,species,max_R,coldens_min=0,coldens_max=1000,minmass=10**11.8,maxmass=10**12.2,savename=None,coldens_vs_R=False,coverfrac_vs_R=False,coverfrac_within_R=False,rudie_LLS=False,rudie_DLA=False):
        # Calculates covering fraction in similar way to observers.  Gather all halos within specified mass range, and treat 
        # all corresponding sightlines together.  

        # Set up bins for radius:
        n_Rbins = 50
        [Rbins_min,Rbins_max] = AU._bin_setup(0.,max_R,n_Rbins)

        for i in np.arange(self.ndat):
            print "Working on {}".format(self.run_list[i])
            # Check if npz for this data has been saved before:
            # foo = "{}kpc".format(max_R)
            # npz_fname = "fc_{}_s{}_{}_N{}_within{}_12minmass.npz".format(self.run_list[i],self.snapnum_list[i],species,Nmin,foo)
            # print self.npz_base+npz_fname
            if False: #os.path.isfile(self.npz_base+npz_fname):
                print "Loaded from npz!"
                f = np.load(self.npz_base+npz_fname)
                m = f['m']
                fc_vs_R = f['fc_vs_R']
                f.close()
            else:
                pass 
                grp_ids = self.find_desired_grpids_withoutsnap(i,minmass=minmass,maxmass=maxmass)
                #gal_dict = self.get_gal_props(i,grp_ids)

                r = np.zeros(0)
                N = np.zeros(0)
                for j in np.arange(np.size(grp_ids)):
                    grp_id = grp_ids[j]
                    grp_data = self.load_grid_data(i,grp_id)
                    grid = grp_data[species]
                    grid_rad = grp_data['grid_radius_pkpc']
                    ngrid = grp_data['ngrid']

                    [gridx,gridy] = np.meshgrid(np.arange(ngrid),np.arange(ngrid))
                    grid_cent = (ngrid-1)/2. #assume square grid: grid_centx = grid_centy = grid_cent
                    r_grid = np.sqrt((gridx-grid_cent)**2+(gridy-grid_cent)**2)
                    r_kpc = self._grid_to_kpc(r_grid,ngrid,grid_rad)
                    
                    dist_thresh = max_R
                    if dist_thresh > grid_rad:
                        raise Exception('Desired radial threshold ({}) exceeded grid radius ({}) for covering fraction calculation.'.format(dist_thresh,grid_rad))

                    r = np.append(r,r_kpc)
                    N = np.append(N,grid)

                if coldens_vs_R:
                    # Calculates the median column density as a function of radius
                    [Q1,med,Q3] = AU._calc_percentiles_v2(r,Rbins_min,Rbins_max,N)

                    plt.close('all')
                    plt.figure()
                    plt.plot(Rbins_min,med,color=self.color_list[i],label=self.label_list)
                    plt.fill_between(Rbins_min,Q1,Q3,color=self.color_list[i],alpha=0.3)
                    plt.xlim([0.,max_R])
                    plt.xlabel('Radius [pkpc]')
                    plt.ylabel(r"log$_{10}$ N$_\mathrm{"+species+"}$ (cm$^{-2}$)")
                    plt.legend()
                    if savename != None:
                        plt.savefig(self.fig_base+savename+".pdf")
                    else:
                        plt.savefig(self.fig_base+"cd_v_R.pdf")

                if coverfrac_within_R:
                    fc = np.zeros(n_Rbins)
                    for k in np.arange(n_Rbins):
                        in_Rbin = r<Rbins_max[k]
                        covered = np.logical_and(N[in_Rbin]>=coldens_min,N[in_Rbin]<coldens_max)
                        fc[k] = np.float(np.sum(covered))/np.float(np.sum(in_Rbin))

                    plt.close('all')
                    plt.figure()
                    plt.plot(Rbins_min,fc,color=self.color_list[i],label=self.label_list)
                    plt.xlim([0.,max_R])
                    plt.xlabel('Radius [pkpc]')
                    plt.ylim([-0.05,1.1])
                    plt.ylabel(r"$f_C (< R)$")

                    # Rudie et al. 2012, cumulative covering fraction of LLS and DLA for M_halo ~ 10^12
                    if rudie_LLS:
                        xdat = np.array([90.,180.])
                        ydat = np.array([0.3,0.24])
                        yerr = np.array([0.14,0.09])
                        xerr = np.array([8.,16.])
                        plt.ylim([0.,1.0])
                        plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,label='Rudie+: LLS',fmt='.',color='purple')
                    elif rudie_DLA:
                        xdat = np.array([90.,180.])
                        ydat = np.array([0.,0.04])
                        yerr = np.array([0.1,0.04])
                        xerr = np.array([8.,16.])
                        plt.ylim([0.,1.0])
                        plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,label='Rudie+: DLA',fmt='.',color='purple')
                    plt.legend()
                    if savename != None:
                        plt.savefig(self.fig_base+savename+".pdf")
                    else:
                        plt.savefig(self.fig_base+"fc_within_R.pdf")


                if coverfrac_vs_R:
                    fc = np.zeros(n_Rbins)
                    for k in np.arange(n_Rbins):
                        in_Rbin = np.logical_and(r<Rbins_max[k],r>Rbins_min[k])
                        covered = np.logical_and(N[in_Rbin]>=coldens_min,N[in_Rbin]<coldens_max)
                        fc[k] = np.float(np.sum(covered))/np.float(np.sum(in_Rbin))

                    plt.close('all')
                    plt.figure()
                    plt.plot(Rbins_min,fc,color=self.color_list[i],label=self.label_list)
                    plt.xlim([0.,max_R])
                    plt.xlabel('Radius [pkpc]')
                    plt.ylim([-0.05,1.1])
                    plt.ylabel(r"$f_C (R)$")

                    plt.legend()
                    if savename != None:
                        plt.savefig(self.fig_base+savename+".pdf")
                    else:
                        plt.savefig(self.fig_base+"fc_vs_R.pdf")








    #############
    # Utilities #
    #############
    def find_desired_grpids_withoutsnap(self,i,minmass=0,maxmass=1e20):
        # usually - read subfind, find desired grp masses, etc.
        # BUT, looks like Mark deleted his snapshots.  In this case, we need to directly use the saved data.

        # get list of snapshot saved data
        grp_ids = np.array([])
        for fn in glob.glob(self.CGMsnap_base+"{}/s{}/*.hdf5".format(self.run_list[i],self.snapnum_list[i])):
            f = h5py.File(fn,'r')
            grp_id = f['Header'].attrs['grp_id']
            grp_mass = AU.PhysicalMass(f['Header'].attrs['grp_mass'])
            f.close()

            if np.logical_and(grp_mass > minmass,grp_mass < maxmass):
                grp_ids = np.append(grp_ids,grp_id)
        return grp_ids

    def load_grid_data(self,i,grp_id):
        # assume that i is singular
        grp_data = {}
        # print self.grid_base+"{}/s{}/{}.hdf5".format(self.run_list[i],self.snapnum_list[i],str(int(grp_id)).zfill(5))
        f = h5py.File(self.grid_base+"{}/s{}/{}.hdf5".format(self.run_list[i],self.snapnum_list[i],str(int(grp_id)).zfill(5)),'r')
        for key in list(f['Header'].attrs): grp_data[key] = f['Header'].attrs[key]
        for key in f['grids'].keys(): grp_data[key] = np.array(f['grids'][key])
        return grp_data 

    def _grid_to_kpc(self,r_grid,ngrid,grid_rad):
        return r_grid * (2.*grid_rad)/(ngrid)


    # self.get_gal_props

    def dat_prep(self,redshifts,c0_128=0,c0_256=0,c0_512=0,c0_fw_256=0,c0_sw_256=0,c2_256=0,c3_512=0,c4_512=0,c4_check=0,g0_BH=0,g10_BH=0,g20_BH=0,g30_BH=0,g40_BH=0,g50_BH=0,g25_BH=0,g75_BH=0,g95_BH=0,g0_noBH=0,g10_noBH=0,g20_noBH=0,g30_noBH=0,g40_noBH=0,g50_noBH=0,g25_noBH=0,g75_noBH=0,g95_noBH=0,g10_nothermal=0,g20_nothermal=0,g30_nothermal=0,g40_nothermal=0,g50_nothermal=0,g25_fixv=0,g50_fixv=0,g25_fixv_fixeta=0,g50_fixv_fixeta=0,g25_fixv_nothermal=0,g50_fixv_nothermal=0,c0_nometalwinds=0,c0_fullmetalwinds=0):
        self.mv_base = "/n/hernquistfs1/mvogelsberger/projects/GFM/Production/Cosmo/"
        self.sb_base = "/n/hernquistfs1/spb/Cosmo/"
        self.js_base = "/n/hernquistfs1/jsuresh/Runs/"

        if c0_256 == 1:
            for redshift in redshifts:
                self.run_list.append('c0_256')
                self.color_list.append('blue')
                self.label_list.append('Fiducial')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(self.mv_base+"Cosmo0_V6/L25_n256")
                if redshift == '4': self.snapnum_list.append(54)
                if redshift == '3': self.snapnum_list.append(60)
                if redshift == '2': self.snapnum_list.append(68)
        if c0_fw_256 == 1:
            for redshift in redshifts:
                self.run_list.append('c0_fw_256')
                self.color_list.append('blue')
                self.label_list.append('Fast Winds')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(self.mv_base+"Cosmo0_V6/L25_n256")
                if redshift == '4': self.snapnum_list.append(54)
                if redshift == '3': self.snapnum_list.append(60)
                if redshift == '2': self.snapnum_list.append(68)
        if c0_sw_256 == 1:
            for redshift in redshifts:
                self.run_list.append('c0_sw_256')
                self.color_list.append('green')
                self.label_list.append('Higher Mass-Loading')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(self.mv_base+"Cosmo0_V6/L25_n256")
                if redshift == '4': self.snapnum_list.append(54)
                if redshift == '3': self.snapnum_list.append(60)
                if redshift == '2': self.snapnum_list.append(68)
        if c0_fw_256 == 1:
            for redshift in redshifts:
                self.run_list.append('c0_fw_256')
                self.color_list.append('blue')
                self.label_list.append('Fast Winds')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(self.mv_base+"Cosmo0_V6/L25_n256")
                if redshift == '4': self.snapnum_list.append(54)
                if redshift == '3': self.snapnum_list.append(60)
                if redshift == '2': self.snapnum_list.append(68)
        if c2_256 == 1:
            for redshift in redshifts:
                self.run_list.append('c2_256')
                self.color_list.append('blue')
                self.label_list.append('No AGN')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(self.mv_base+"Cosmo2_V6/L25_n256")
                if redshift == '4': self.snapnum_list.append(54)
                if redshift == '3': self.snapnum_list.append(60)
                if redshift == '2': self.snapnum_list.append(68)
        if c2_256 == 1:
            for redshift in redshifts:
                self.run_list.append('c2_256')
                self.color_list.append('cyan')
                self.label_list.append('No AGN')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(self.mv_base+"Cosmo2_V6/L25_n256")
                if redshift == '4': self.snapnum_list.append(54)
                if redshift == '3': self.snapnum_list.append(60)
                if redshift == '2': self.snapnum_list.append(68)

        if g0_BH == 1:
            for redshift in redshifts:
                self.run_list.append('c0_256')
                self.color_list.append('blue')
                self.label_list.append('Fiducial')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(self.mv_base+"Cosmo0_V6/L25_n256")
                if redshift == '4': self.snapnum_list.append(54)
                if redshift == '3': self.snapnum_list.append(60)
                if redshift == '2': self.snapnum_list.append(68)

        if g50_fixv_nothermal == 1:
            for redshift in redshifts:
                self.run_list.append('g50_fixv_nothermal')
                self.color_list.append('chartreuse')
                self.label_list.append('Faster Winds')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(self.js_base+"gam_50_fixv_nothermal/")
                if redshift == '4': self.snapnum_list.append(1)
                if redshift == '3': self.snapnum_list.append(3)
                if redshift == '2': self.snapnum_list.append(5)
        if g50_BH == 1:
            for redshift in redshifts:
                self.run_list.append('g50_BH')
                self.color_list.append('brown')
                self.label_list.append('Fixed-E Hot Winds')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(self.js_base+"gam_50_BH/")
                if redshift == '4': self.snapnum_list.append(1)
                if redshift == '3': self.snapnum_list.append(3)
                if redshift == '2': self.snapnum_list.append(5)
        if g50_fixv == 1:
            for redshift in redshifts:
                self.run_list.append('g50_fixv')
                self.color_list.append('magenta')
                self.label_list.append('Fixed-v Hot Winds')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(self.js_base+"gam_50_fixv/")
                if redshift == '4': self.snapnum_list.append(1)
                if redshift == '3': self.snapnum_list.append(3)
                if redshift == '2': self.snapnum_list.append(5)
        if c0_nometalwinds == 1:
            for redshift in redshifts:
                self.run_list.append('c0_nometalwinds')
                self.color_list.append('gray')
                self.label_list.append('Pristine Winds')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(self.js_base+"c0_nometalwinds/")
                if redshift == '4': self.snapnum_list.append(1)
                if redshift == '3': self.snapnum_list.append(3)
                if redshift == '2': self.snapnum_list.append(5)
        if c0_fullmetalwinds == 1:
            for redshift in redshifts:
                self.run_list.append('c0_fullmetalwinds')
                self.color_list.append('black')
                self.label_list.append('Fully Enriched Winds')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(self.js_base+"c0_fullmetalwinds/")
                if redshift == '4': self.snapnum_list.append(1)
                if redshift == '3': self.snapnum_list.append(3)
                if redshift == '2': self.snapnum_list.append(5)

        self.ndat = np.size(self.run_list)

if __name__ == '__main__':
    full_analysis()

