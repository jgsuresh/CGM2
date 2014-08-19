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

import SphMap as sm 


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



        # # Cosmo0_V6 (512^3 particles):
        # redshifts = ['2']
        # self.dat_prep(redshifts,c0_256=1)
        # redshifts = ['2']
        # self.dat_prep(redshifts,c0_512=1)
        
        # # Cosmo3_V6 (512^3 particles):
        # redshifts = ['4','3','2']
        # self.dat_prep(redshifts,c3_512=1)

        # # Cosmo4_V6 (512^3 particles):
        # redshifts = ['2.5']
        # self.dat_prep(redshifts,c4_512=1)
        # self.dat_prep(redshifts,c4_check=1)

        # # # Cosmo2_V6 (256^3 particles):
        # redshifts = ['2']
        # self.dat_prep(redshifts,c2_256=1)

        # # Cosmo0_V6_fastWinds (256^3 particles):
        # redshifts = ['4','3','2','1','0.3','0']
        # self.dat_prep(redshifts,c0_fw_256=1)


        # # Cosmo5_V6 (256^3 particles):
        # redshifts = ['4','3','2']
        # self.dat_prep(redshifts,c5_256=1)

        # Gamma boxes (256^3 particles):
        # redshifts = ['2']
        # self.dat_prep(redshifts,g0_BH=1)
        # # Cosmo2_V6 (256^3 fparticles):
        # redshifts = ['2']
        # self.dat_prep(redshifts,c2_256=1)

        # # # # Cosmo0_V6_strongWinds (256^3 particles):
        # redshifts = ['2']
        # self.dat_prep(redshifts,c0_sw_256=1)

        # redshifts = ['2']
        # self.dat_prep(redshifts,g50_fixv_nothermal=1)

        # redshifts = ['2']
        # self.dat_prep(redshifts,g10_BH=1)
        # redshifts = ['2']
        # self.dat_prep(redshifts,g20_BH=1)
        # redshifts = ['2']
        # self.dat_prep(redshifts,g30_BH=1)
        # redshifts = ['2']
        # self.dat_prep(redshifts,g40_BH=1)
        # redshifts = ['2']
        # self.dat_prep(redshifts,g50_BH=1)
        # redshifts = ['0.3']
        # # self.dat_prep(redshifts,g25_BH=1)
        # # redshifts = ['2']
        # # self.dat_prep(redshifts,g25_noBH=1)
        # # redshifts = ['3']
        # # self.dat_prep(redshifts,g75_BH=1)
        # # redshifts = ['4','3','2']
        # # self.dat_prep(redshifts,g75_noBH=1)
        # # redshifts = ['3']
        # # self.dat_prep(redshifts,g95_BH=1)
        # # redshifts = ['4','3','2']
        # # self.dat_prep(redshifts,g95_noBH=1)

        # redshifts = ['2']
        # self.dat_prep(redshifts,g10_noBH=1)
        # self.dat_prep(redshifts,g20_noBH=1)
        # self.dat_prep(redshifts,g30_noBH=1)
        # self.dat_prep(redshifts,g40_noBH=1)
        # self.dat_prep(redshifts,g50_noBH=1)

        # redshifts = ['2']
        # self.dat_prep(redshifts,g10_nothermal=1)
        # redshifts = ['2']
        # self.dat_prep(redshifts,g20_nothermal=1)
        # redshifts = ['2']
        # self.dat_prep(redshifts,g30_nothermal=1)
        # redshifts = ['2']
        # self.dat_prep(redshifts,g40_nothermal=1)
        # redshifts = ['2']
        # self.dat_prep(redshifts,g50_nothermal=1)


        # redshifts = ['2']
        # self.dat_prep(redshifts,g25_fixv=1)
        # redshifts = ['2']
        # self.dat_prep(redshifts,g50_fixv=1)

        # redshifts = ['2']
        # self.dat_prep(redshifts,g25_fixv_fixeta=1)
        # redshifts = ['2']
        # self.dat_prep(redshifts,g50_fixv_fixeta=1)

        # redshifts = ['2']
        # self.dat_prep(redshifts,g25_fixv_nothermal=1)
        # redshifts = ['2']
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
        # self.enrichment_vs_volume()

        ##############################
        # 3D data analysis functions #
        ##############################
        # self.los_profile("metal",mmin=10.**11.,mmax=10**11.5,profile_type=2,savename="avg_metal_z03_fixv")
        # self.los_profile("T",mmin=10.**11.,mmax=10**11.5,profile_type=2,savename="avg_T_z03_fixv")
        # self.los_profile("O",mmin=10.**11.,mmax=10**11.5,profile_type=2,savename="O_z03_fixv")
        # self.galaxy_mass_hist()
        # self.phase_diagram(mmin=10.**11,savename="phase_c0")
        # self.phase_diagram_diff(savename="phasediff_c0_check")
        # self.phase_budget('m',savename='m_budget_lowmetal',minimal=False)
        # self.phase_budget('z',savename='z_budget_lowmetal',minimal=False)
        # self.plot_individual_halos(grp_ids=np.array([24, 25, 28, 29, 35, 41, 46, 63, 64]),true_hsml=True) #c0_256
        # self.plot_individual_halos(grp_ids=np.array([22, 28, 25, 26, 31, 52, 45, 61, 76]),true_hsml=True) #gam_50_BH
        # self.plot_individual_halos(grp_ids=np.array([23, 36, 28, 29, 37, 53, 48, 71, 68]),true_hsml=True) #gam_50_fixv
        # self.plot_individual_halos(grp_ids=np.array([23, 30, 26, 28, 35, 55, 50, 71, 63]),true_hsml=True) #gam_50_fixv_nothermal
        # self.CGM_and_gal_metallicity(savename="massmet_test")
        # self.radial_profile('T')
        # self.radial_profile('z')

        ##############################
        # 2D data analysis functions #
        ##############################
        # self.cover_frac_within_R("H1",17,kpc_thresh=180.,rud_pro=True,plot_mode=3,savename="HI_fc_z2_90kpc_binned")
        # self.cover_frac_within_R("O6",14,kpc_thresh=100.,plot_mode=3,savename="O6_fc_z3_100kpc_binned")
        # self.cover_frac_vs_R("H1",17.2,200.,Nmax=20.3,minmass=10**11.8,maxmass=10**12.2,plot_mode=1,savename="LLS_fiduc",rudie_LLS=True)
        # self.cover_frac_vs_R("H1",20.3,200.,minmass=10**11.8,maxmass=10**12.2,plot_mode=1,savename="DLA_c0",rudie_DLA=True)
        # self.observational_cover_frac_vs_R("H1",17.2,200.,Nmax=20.3,minmass=10**11.8,maxmass=10**12.2,savename="LLS_postisaac",rudie_LLS=True)
        # self.observational_cover_frac_vs_R("H1",20.3,200.,minmass=10**11.8,maxmass=10**12.2,savename="DLA_postisaac",rudie_DLA=True)
        # self.cover_frac_vs_R("Mg2",12.5,200.,minmass=10**12.0,plot_mode=1,savename="MgII_z03_gam_12")
        # self.cover_frac_vs_R("Si3",12.75,200.,minmass=10**12.0,plot_mode=1,savename="SiIII_z03_gam_12")
        # self.cover_frac_vs_R("O6",14,200.,minmass=10**11,maxmass=10**15,plot_mode=1,savename="O6_z2_gam_11")
        # self.columndens_vs_R("O6",150.,minmass=10**11.8,maxmass=10**12.2,plot_mode=1,savename="O6_z3_150kpc_posti")
        # self.plot_grids("H1",minmass=10.**11.8,vmax=25.)
        # self.plot_grids("O6",minmass=10.**11.5,vmax=16.)
        # self.COS_Halos_comparison("Mg2",savename="Mg2_Werk_fixv")
        # self.COS_Halos_comparison("Si3",savename="Si3_Werk_fixv")
        # self.COS_Halos_comparison("O6",savename="O6_Werk_fixv")

        # OTHER:
        # print "Virial Temperature: ",AU.HaloTvir(10.**12.,2.)

        

    ############################################################################################################
    ############################################################################################################
    
    def cover_frac_within_R(self,species,N_thresh,Rvir_thresh=-1,kpc_thresh=-1,rud_pro=False,plot_mode=1,run_mode='save',savename=None):
        # Plots covering fraction within a set distance R.  

        #Plot mode is:
        # 1: Simple scatter plot of covering fraction vs group mass
        # 2: Mass binned, with only median line shown
        # 3: Mass binned, with median and scatter about median shown

        plt.close('all')
        plt.figure()

        if Rvir_thresh > -1 and kpc_thresh == -1: use_Rvir = True
        if Rvir_thresh == -1 and kpc_thresh > -1: use_Rvir = False
        if Rvir_thresh > -1 and kpc_thresh > -1:
            print "Can't specify both Rvir and kpc for covering fraction!"

        for i in np.arange(self.ndat):
            # Check if npz for this data has been saved before:
            if use_Rvir: foo = "{}Rv".format(int(Rvir_thresh))
            else: foo = "{}kpc".format(kpc_thresh)
            npz_fname = "fc_{}_s{}_{}_N{}_{}.npz".format(self.run_list[i],self.snapnum_list[i],species,N_thresh,foo)
            print self.npz_base+npz_fname
            if os.path.isfile(self.npz_base+npz_fname):
                print "Loaded from npz!"
                f = np.load(self.npz_base+npz_fname)
                m = f['m']
                fc = f['fc']
                f.close()
            else:
                grp_ids = self.find_desired_grpids(i,minmass=1e11)
                m = np.zeros(np.size(grp_ids))
                fc = np.zeros(np.size(grp_ids))
                for grp_id in grp_ids:
                    grp_data = self.load_grid_data(i,grp_id)
                    grid = grp_data[species]
                    grid_rad = grp_data['grid_radius_pkpc']
                    ngrid = grp_data['ngrid']

                    [gridx,gridy] = np.meshgrid(np.arange(ngrid),np.arange(ngrid))
                    grid_cent = (ngrid-1)/2. #assume square grid: grid_centx = grid_centy = grid_cent
                    r_grid = np.sqrt((gridx-grid_cent)**2+(gridy-grid_cent)**2)
                    r_kpc = self._grid_to_kpc(r_grid,ngrid,grid_rad)

                    if use_Rvir:
                        grp_Rvir = grp_data['grp_Rvir']
                        dist_thresh = Rvir_thresh*grp_Rvir
                    else:
                        dist_thresh = kpc_thresh
                    if dist_thresh > grid_rad:
                        raise Exception('Desired radial threshold ({}) exceeded grid radius ({}) for covering fraction calculation.'.format(dist_thresh,grid_rad))
                    within = (r_kpc <= dist_thresh)

                    n_tot = np.sum(within)
                    n_cov = np.sum(grid[within] > N_thresh)

                    m = np.append(m,AU.PhysicalMass(grp_data['grp_mass']))
                    fc = np.append(fc,np.float(n_cov)/np.float(n_tot))

                # Save covering fraction data to npz file
                np.savez(self.npz_base+npz_fname,run=self.run_list[i],snapnum=self.snapnum_list[i],species=species,N_thresh=N_thresh,Rvir_thresh=Rvir_thresh,kpc_thresh=kpc_thresh,grp_ids=grp_ids,m=m,fc=fc)  


            if plot_mode == 1:
                plt.scatter(np.log10(m),fc,color=self.color_list[i],label=self.label_list[i])
                plt.legend()

            elif plot_mode == 2 or plot_mode == 3:
                x = np.log10(m)
                [binmin,binmax] = AU._bin_setup(11.,np.max(x),10,logbins=True)
                binmin = np.log10(binmin)
                binmax = np.log10(binmax)
                binmed = (binmin+binmax)/2.
                [Q1,med,Q3] = AU._calc_percentiles_v2(x,binmin,binmax,fc,min_percentile=10,max_percentile=90)

                plt.plot(binmed,med,color=self.color_list[i],label=self.label_list[i],linestyle=self.linestyle_list[i])
                plt.legend()
                if plot_mode == 3:
                    plt.fill_between(binmed,Q3,Q1,color=self.color_list[i],alpha=0.3)


        if rud_pro:
            if use_Rvir:
                if Rvir_thresh == 1.:
                    plt.errorbar(np.array([12.]),np.array([0.3]),xerr=0.2,yerr=0.19,marker='s',color='purple',label='Rudie+')
                    plt.errorbar(np.array([12.5]),np.array([0.65]),xerr=0.2,yerr=0.10,marker='^',color='purple',label='Prochaska+')
                    plt.ylim([-0.05,1.1])
                    plt.ylabel(r"$f_C(< 1 R_{200})$")
                    plt.legend(prop={'size':6})
                elif Rvir_thresh == 2.:
                    plt.errorbar(np.array([12.]),np.array([0.32]),xerr=0.2,yerr=0.13,marker='s',color='purple',label='Rudie+')
                    plt.errorbar(np.array([12.5]),np.array([0.45]),xerr=0.2,yerr=0.10,marker='^',color='purple',label='Prochaska+')
                    plt.ylim([-0.02,0.6])
                    plt.ylabel(r"$f_C(< 2 R_{200})$")
                    plt.legend(prop={'size':6})

            else:
                if kpc_thresh == 90:
                    plt.errorbar(np.array([12.]),np.array([0.3]),xerr=0.2,yerr=0.19,marker='s',color='purple',label='Rudie+')
                    plt.errorbar(np.array([12.5]),np.array([0.65]),xerr=0.15,yerr=0.10,marker='^',color='purple',label='Prochaska+')
                    plt.ylim([-0.05,1.1])
                    plt.ylabel(r"$f_C(< 90$ kpc$)$")
                    plt.legend(prop={'size':6})
                elif kpc_thresh == 180:
                    plt.errorbar(np.array([12.]),np.array([0.32]),xerr=0.2,yerr=0.13,marker='s',color='purple',label='Rudie+')
                    plt.errorbar(np.array([12.5]),np.array([0.45]),xerr=0.2,yerr=0.10,marker='^',color='purple',label='Prochaska+')
                    plt.ylim([-0.02,0.6])
                    plt.ylabel(r"$f_C(< 180$ kpc$)$")
                    plt.legend(prop={'size':6})

        if savename != None:
            plt.savefig(self.fig_base+savename+".pdf")
        else:
            plt.savefig(self.fig_base+"fc_within_R.pdf")




    def cover_frac_vs_R(self,species,Nmin,kpc_thresh,Nmax=1000.,minmass=0.,maxmass=1e18,rud_pro=False,plot_mode=1,run_mode='save',savename=None,rudie_LLS=False,rudie_DLA=False):
        # Plots cumulative radial covering fraction for a given range of halo masses.  

        #Plot mode is:
        # 1: Only median line shown
        # 2: Median and scatter about median shown

        plt.close('all')
        plt.figure()

        # Set up bins for radius:
        n_Rbins = 50
        [Rbins_min,Rbins_max] = AU._bin_setup(0.,200.,n_Rbins)

        for i in np.arange(self.ndat):
            # Check if npz for this data has been saved before:
            foo = "{}kpc".format(kpc_thresh)
            npz_fname = "fc_{}_s{}_{}_N{}_within{}_12minmass.npz".format(self.run_list[i],self.snapnum_list[i],species,Nmin,foo)
            print self.npz_base+npz_fname
            if os.path.isfile(self.npz_base+npz_fname):
                print "Loaded from npz!"
                f = np.load(self.npz_base+npz_fname)
                m = f['m']
                fc_vs_R = f['fc_vs_R']
                f.close()
            else:
                grp_ids = self.find_desired_grpids(i,minmass=minmass,maxmass=maxmass)
                #gal_dict = self.get_gal_props(i,grp_ids)

                print "grp_ids: ",grp_ids
                print "np.size(grp_ids) ",np.size(grp_ids)
                m = np.zeros(np.size(grp_ids))
                fc = np.zeros(np.size(grp_ids))
                fc_vs_R = np.zeros([np.size(grp_ids),n_Rbins])
                for j in np.arange(np.size(grp_ids)):
                    grp_id = grp_ids[j]
      