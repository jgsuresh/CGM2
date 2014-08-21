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
        # redshifts = ['2']
        # self.dat_prep(redshifts,g0_BH=1)
        # self.dat_prep(redshifts,c2_256=1)
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
        print self.snapdir_list
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
        self.CGM_and_gal_metallicity(savename="massmet_newfan")
        # self.radial_profile('T')
        # self.radial_profile('z')

        ##############################
        # 2D data analysis functions #
        ##############################
        # self.grid_sightlines("H1",200.,coldens_min=15.5,minmass=10**11.8,maxmass=10**12.2,coverfrac_within_R=True,rudie_155=True,savename="rudie_155_test",show_Fumagalli=True)
        # self.grid_sightlines("H1",200.,coldens_min=17.2,minmass=10**11.8,maxmass=10**12.2,coverfrac_within_R=True,rudie_172=True,savename="rudie_172_z2",show_Fumagalli=True)
        # self.grid_sightlines("H1",200.,coldens_min=19.,minmass=10**11.8,maxmass=10**12.2,coverfrac_within_R=True,rudie_19=True,savename="rudie_19_z2",show_Fumagalli=True)
        # self.grid_sightlines("H1",200.,coldens_min=20.3,minmass=10**11.8,maxmass=10**12.2,coverfrac_within_R=True,rudie_203=True,savename="rudie_203_z2",show_Fumagalli=True)
        # print "2"
        # self.grid_sightlines("H1",200.,coldens_min=20.3,minmass=10**11.8,maxmass=10**12.2,coverfrac_within_R=True,rudie_DLA=True,savename="DLA_newfan2")
        # print "3"
        # self.grid_sightlines("H1",200.,minmass=10**11.8,maxmass=10**12.2,coldens_vs_R=True,savename="coldens_newfan2")

    ############################################################################################################
    ############################################################################################################

    #########################################
    # Full-snapshot data analysis functions #
    #########################################

    ##############################
    # 3D data analysis functions #
    ##############################

    def CGM_and_gal_metallicity(self,savename=None,plot_mode=1):
        # Compute the mass and metal content of the CGM:

        plt.close('all')    
        plt.figure(figsize=(3.54331*2,3.14))

        p1 = plt.subplot(1,2,1)
        p2 = plt.subplot(1,2,2)
        subplot_list = [p1,p2]

        lines = []

        # Bin halos by mass
        mbins_min = 10.**(10.+np.array([1.0,1.1,1.2,1.3,1.4,1.5,1.7,1.9,2.1,2.4]))
        mbins_max = 10.**(10.+np.array([1.1,1.2,1.3,1.4,1.5,1.7,1.9,2.1,2.4,2.7]))
        mbins_med = 10.**((np.log10(mbins_min) + np.log10(mbins_max))/2.)
        n_mbins = np.size(mbins_min)

        for i in np.arange(self.ndat):
            print "Working on {}".format(self.run_list[i]) 

            npz_fname = "CGMmassmet_{}_s{}.npz".format(self.run_list[i],self.snapnum_list[i])

            if os.path.isfile(self.npz_base+npz_fname):
                print "Loaded from npz!: {}".format(self.npz_base+npz_fname)
                f = np.load(self.npz_base+npz_fname)
                grp_m = f['grp_m']
                full_dat = f['full_dat']
                f.close()
            else: 
                cat = readsubfHDF5.subfind_catalog(self.snapdir_list[i],self.snapnum_list[i],long_ids=True)
                grp_ids = self.find_desired_grpids(i,minmass=1e11)

                m_CGM = np.zeros(0)
                m_ISM = np.zeros(0)
                m_stars = np.zeros(0)
                zm_CGM = np.zeros(0)
                zm_ISM = np.zeros(0)
                zm_stars = np.zeros(0)
                grp_m = np.zeros(0)
                for j in np.arange(np.size(grp_ids)):
                    grp_id = grp_ids[j]
                    grp_data = self.load_CGMsnap_data(i,grp_id)
                    grp_m = np.append(grp_m,AU.PhysicalMass(grp_data['grp_mass']))
                    gas_m = AU.PhysicalMass(np.array(grp_data['Masses']))
                    z = np.array(grp_data['GFM_Metallicity'])
                    r = AU.PhysicalPosition(self._calc_radii(grp_data),grp_data['redshift'])
                    grp_Rvir = AU.PhysicalPosition(grp_data['grp_Rvir'],grp_data['redshift'])

                    print "(Max particle radius)/(3 Rvir) = {}".format(np.max(r)/(3*grp_Rvir))
                    r_ind = r <= 3.*grp_Rvir

                    # Find ISM gas (i.e. gas that is star-forming)
                    SF = self._find_ISM_gas(grp_data)
                    ISM = np.logical_and(r_ind,SF)
                    CGM = np.logical_and(r_ind,np.logical_not(SF))

                    m_ISM = np.append(m_ISM,np.sum(gas_m[ISM]))
                    m_CGM = np.append(m_CGM,np.sum(gas_m[CGM]))
                    m_stars = np.append(m_stars,AU.PhysicalMass(cat.GroupMassType[grp_id,4]))
                    zm_ISM = np.append(zm_ISM,np.sum(gas_m[ISM]*z[ISM]))
                    zm_CGM = np.append(zm_CGM,np.sum(gas_m[CGM]*z[CGM]))
                    zm_stars = np.append(zm_stars,AU.PhysicalMass(cat.GroupMassType[grp_id,4])*cat.GroupStarMetallicity[grp_id])

                full_dat = [m_ISM+m_stars,m_CGM,zm_ISM+zm_stars,zm_CGM]

                np.savez(self.npz_base+npz_fname,full_dat=full_dat,grp_m=grp_m)


            med = np.zeros([n_mbins,4])
            if plot_mode == 2:
                Q1 = np.zeros([n_mbins,4])
                Q3 = np.zeros([n_mbins,4])

            for mbin_ind in np.arange(n_mbins):
                mmin = mbins_min[mbin_ind]
                mmax = mbins_max[mbin_ind]
                m_select = np.logical_and(grp_m > mmin, grp_m <= mmax)
                for j in np.arange(4):
                    med[mbin_ind,j] = np.median(full_dat[j][m_select])
                    if plot_mode == 2:
                        Q1[mbin_ind,j] = np.percentile(full_dat[j][m_select],10)
                        Q3[mbin_ind,j] = np.percentile(full_dat[j][m_select],90)

            if i==0:
                p_gal,=p1.semilogy(np.log10(mbins_med),med[:,0]/((0.0456/0.27)*mbins_med),color=self.color_list[i],linestyle='dashed')
                p_CGM,=p1.semilogy(np.log10(mbins_med),med[:,1]/((0.0456/0.27)*mbins_med),color=self.color_list[i],label=self.label_list[i])
                legend1 = p1.legend([p_gal,p_CGM], ["Stars+ISM","CGM"], prop={'size':6})
                plt.gca().add_artist(legend1)
            else:
                # p1.plot(np.log10(mbins_med),np.log10(med[:,0]/((0.0456/0.27)*mbins_med)),color=self.color_list[i],linestyle='dotted')
                # p1.plot(np.log10(mbins_med),np.log10(med[:,1]/((0.0456/0.27)*mbins_med)),color=self.color_list[i],label=self.label_list[i])
                p1.semilogy(np.log10(mbins_med),med[:,0]/((0.0456/0.27)*mbins_med),color=self.color_list[i],linestyle='dashed')
                p = p1.semilogy(np.log10(mbins_med),med[:,1]/((0.0456/0.27)*mbins_med),color=self.color_list[i],label=self.label_list[i])
                lines.append(p)

            xdat1 = (med[:,2]/med[:,0])/AU.SolarMetallicity
            xdat2 = (med[:,3]/med[:,1])/AU.SolarMetallicity
            # p2.plot(np.log10(mbins_med),np.log10(xdat1),color=self.color_list[i],linestyle='dotted')
            # p2.plot(np.log10(mbins_med),np.log10(xdat2),color=self.color_list[i],label=self.label_list[i])
            p2.semilogy(np.log10(mbins_med),xdat1,color=self.color_list[i],linestyle='dashed')
            p2.semilogy(np.log10(mbins_med),xdat2,color=self.color_list[i],label=self.label_list[i])

        #p1.set_ylim([-1.,1.])
        p1.set_ylim([0.1,10.])
        p1.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        p1.xaxis.set_ticks(np.arange(11.0,12.6,0.5))
        p1.set_xlabel(r"$\log_{10}\left[M_{\rm Halo}\right]$")
        # p1.set_ylabel(r"$\log_{10}\left[ \frac{M_{\rm baryons}}{M_{\rm Halo} \times \Omega_b} \right]$")
        p1.set_ylabel(r"$\frac{M_{\rm baryons}}{M_{\rm Halo} \times f_b}$")
        p1.legend(loc=2,prop={'size':6})
        # plt.gca().add_artist(legend1)

        p2.xaxis.set_ticks(np.arange(11.0,12.6,0.5))
        # p2.set_ylabel(r"$\log_{10}\left[ \frac{Z_{\rm baryons}}{Z_\odot} \right]$")
        p2.set_ylabel(r"$\frac{Z_{\rm baryons}}{Z_\odot}$")
        p2.set_xlabel(r"$\log_{10}\left[M_{\rm Halo}\right]$")
        #p2.legend(loc=2,prop={'size':5})

        plt.subplots_adjust(wspace=0.4)

        if savename != None: 
            plt.savefig(self.fig_base+savename+".pdf", bbox_inches='tight')
        else:
            plt.savefig(self.fig_base+"mass_budget.pdf", bbox_inches='tight')



    ##############################
    # 2D data analysis functions #
    ##############################

    def grid_sightlines(self,species,max_R,coldens_min=0,coldens_max=1000,minmass=10**11.9,maxmass=10**12.1,savename=None,coldens_vs_R=False,coverfrac_vs_R=False,coverfrac_within_R=False,rudie_155=False,rudie_172=False,rudie_19=False,rudie_203=False,show_Fumagalli=False):
        # Calculates covering fraction in similar way to observers.  Gather all halos within specified mass range, and treat 
        # all corresponding sightlines together.  

        plt.close('all')
        plt.figure(figsize=(3.54331,3.14))

        # Set up bins for radius:
        n_Rbins = 50
        [Rbins_min,Rbins_max] = AU._bin_setup(0.,max_R,n_Rbins)
        Rbins_med = (Rbins_min+Rbins_max)/2.

        for i in np.arange(self.ndat):
            print "Working on {}".format(self.run_list[i])
            
            grp_ids = self.find_desired_grpids(i,minmass=minmass,maxmass=maxmass)
            grp_ids2 = self.find_desired_grpids_withoutsnap(i,minmass=minmass,maxmass=maxmass)
            if not np.array_equal(grp_ids,grp_ids2):
                print "somethings wrong"
            # print "grp_ids ",np.sort(grp_ids)

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

                plt.plot(Rbins_min,med,color=self.color_list[i],label=self.label_list[i])
                plt.fill_between(Rbins_med,Q1,Q3,color=self.color_list[i],alpha=0.3)

            elif coverfrac_vs_R:
                fc = np.zeros(n_Rbins)
                for k in np.arange(n_Rbins):
                    in_Rbin = np.logical_and(r<Rbins_max[k],r>Rbins_min[k])
                    covered = np.logical_and(N[in_Rbin]>=coldens_min,N[in_Rbin]<coldens_max)
                    fc[k] = np.float(np.sum(covered))/np.float(np.sum(in_Rbin))

                plt.plot(Rbins_med,fc,color=self.color_list[i],label=self.label_list[i])
                
            elif coverfrac_within_R:
                fc = np.zeros(n_Rbins)
                for k in np.arange(n_Rbins):
                    in_Rbin = r<Rbins_max[k]
                    covered = np.logical_and(N[in_Rbin]>=coldens_min,N[in_Rbin]<coldens_max)
                    fc[k] = np.float(np.sum(covered))/np.float(np.sum(in_Rbin))

                plt.plot(Rbins_med,fc,color=self.color_list[i],label=self.label_list[i])


        # Rudie et al. 2012, cumulative covering fraction for different values of N_HI for M_halo ~ 10^12
        if rudie_155:
            xdat = np.array([90.,180.])
            ydat = np.array([0.9,0.6])
            yerr = np.array([0.09,0.1])
            xerr = np.array([8.,16.])
            plt.ylim([0.,1.6])
            plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,fmt='.',color='purple')
            if show_Fumagalli:
                plt.scatter(np.array([90.,180.]),np.array([0.38,0.22]),marker='D',color='gray')
        if rudie_172:
            xdat = np.array([90.,180.])
            ydat = np.array([0.3,0.28])
            yerr = np.array([0.14,0.09])
            xerr = np.array([8.,16.])
            plt.ylim([0.,1.0])
            plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,fmt='.',color='purple')
            if show_Fumagalli:
                plt.scatter(np.array([90.,180.]),np.array([0.16,0.07]),marker='D',color='gray')
        elif rudie_19:
            xdat = np.array([90.,180.])
            ydat = np.array([0.1,0.08])
            yerr = np.array([0.09,0.05])
            xerr = np.array([8.,16.])
            plt.ylim([0.,0.4])
            plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,fmt='.',color='purple')
            if show_Fumagalli:
                plt.scatter(np.array([90.,180.]),np.array([0.06,0.03]),marker='D',color='gray')
        elif rudie_203:
            xdat = np.array([90.,180.])
            ydat = np.array([0.,0.04])
            yerr = np.array([0.1,0.04])
            xerr = np.array([8.,16.])
            plt.ylim([0.,0.4])
            plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,fmt='.',color='purple')
            if show_Fumagalli:
                plt.scatter(np.array([90.,180.]),np.array([0.03,0.01]),marker='D',color='gray')


            
        if coldens_vs_R:
            plt.xlim([0.,max_R])
            plt.xlabel('Radius [pkpc]')
            plt.ylabel(r"log$_{10}$ N$_\mathrm{"+species+"}$ (cm$^{-2}$)")
            plt.legend(prop={'size':5},ncol=2)
            plt.subplots_adjust(left=0.2,bottom=0.18)
            if savename != None:
                plt.savefig(self.fig_base+savename+".pdf", bbox_inches='tight')
            else:
                plt.savefig(self.fig_base+"cd_v_R.pdf", bbox_inches='tight')
        elif coverfrac_vs_R:
            plt.xlim([0.,max_R])
            plt.xlabel('Radius [pkpc]')
            plt.ylim([-0.05,1.1])
            plt.ylabel(r"$f_C (R)$")
            plt.legend(prop={'size':5},ncol=2)
            plt.subplots_adjust(left=0.2,bottom=0.18)
            if savename != None:
                plt.savefig(self.fig_base+savename+".pdf", bbox_inches='tight')
            else:
                plt.savefig(self.fig_base+"fc_vs_R.pdf", bbox_inches='tight')
        elif coverfrac_within_R:
            plt.xlim([0.,max_R])
            plt.xlabel('Radius [pkpc]')
            plt.ylabel(r"$f_C (< R)$")
            plt.legend(prop={'size':5},ncol=2)
            plt.subplots_adjust(left=0.2,bottom=0.18)
            if savename != None:
                plt.savefig(self.fig_base+savename+".pdf", bbox_inches='tight')
            else:
                plt.savefig(self.fig_base+"fc_within_R.pdf", bbox_inches='tight')

                



    #############
    # Utilities #
    #############
    def find_desired_grpids_withoutsnap(self,i,minmass=0,maxmass=1e20):
        # usually - read subfind, find desired grp masses, etc.
        # BUT, looks like Mark deleted his snapshots.  In this case, we need to directly use the saved data.

        # get list of snapshot saved data
        grp_ids = np.array([])
        grp_masses = np.array([])
        for fn in glob.glob(self.CGMsnap_base+"{}/s{}/*.hdf5".format(self.run_list[i],self.snapnum_list[i])):
            f = h5py.File(fn,'r')
            grp_id = f['Header'].attrs['grp_id']
            grp_mass = AU.PhysicalMass(f['Header'].attrs['grp_mass'])
            f.close()

            if np.logical_and(grp_mass > minmass,grp_mass < maxmass):
                grp_ids = np.append(grp_ids,grp_id)
                grp_masses = np.append(grp_masses,grp_mass)
        # for j in np.arange(np.size(grp_ids)):
        #     print "grp id {} - grp mass {}".format(grp_ids[j],grp_masses[j])

        return np.sort(grp_ids)

    def find_desired_grpids(self,i,minmass=0,maxmass=1e20):
        print self.snapdir_list[i],self.snapnum_list[i]
        cat = readsubfHDF5.subfind_catalog(self.snapdir_list[i],self.snapnum_list[i],subcat=False)
        grp_mass = AU.PhysicalMass(np.array(cat.Group_M_Crit200))
        mass_select = np.logical_and(grp_mass > minmass,grp_mass < maxmass)

        grp_ids = np.arange(np.float64(cat.ngroups))
        grp_ids = grp_ids[mass_select]
        # grp_mass = grp_mass[mass_select]
        # for j in np.arange(np.size(grp_ids)):
        #     print "grp id {} - grp mass {}".format(grp_ids[j],grp_mass[j])
        return np.sort(grp_ids)

    def load_grid_data(self,i,grp_id):
        # assume that i is singular
        grp_data = {}
        # print self.grid_base+"{}/s{}/{}.hdf5".format(self.run_list[i],self.snapnum_list[i],str(int(grp_id)).zfill(5))
        try: f = h5py.File(self.grid_base+"{}/s{}/{}.hdf5".format(self.run_list[i],self.snapnum_list[i],str(int(grp_id)).zfill(5)),'r')
        except: 
            print "File could not be opened!: {}".format(self.grid_base+"{}/s{}/{}.hdf5".format(self.run_list[i],self.snapnum_list[i],str(int(grp_id)).zfill(5)))
        for key in list(f['Header'].attrs): grp_data[key] = f['Header'].attrs[key]
        for key in f['grids'].keys(): grp_data[key] = np.array(f['grids'][key])
        return grp_data 

    def _grid_to_kpc(self,r_grid,ngrid,grid_rad):
        return r_grid * (2.*grid_rad)/(ngrid)


        
    # self.get_gal_props


    def dat_prep(self,redshifts,c0_256=0,c0_512=0,c3_512=0,c4_512=0,c4_check=0,c2_256=0,c0_fw_256=0,c0_sw_256=0,c5_256=0,g0_BH=0,g10_BH=0,g20_BH=0,g25_BH=0,g30_BH=0,g40_BH=0,g25_noBH=0,g50_BH=0,g75_BH=0,g75_noBH=0,g95_BH=0,g95_noBH=0,g10_noBH=0,g20_noBH=0,g30_noBH=0,g40_noBH=0,g50_noBH=0,g10_nothermal=0,g20_nothermal=0,g30_nothermal=0,g40_nothermal=0,g50_nothermal=0,g25_fixv=0,g50_fixv=0,g25_fixv_nothermal=0,g50_fixv_nothermal=0,g25_fixv_fixeta=0,g50_fixv_fixeta=0,c0_nometalwinds=0,c0_fullmetalwinds=0):
        mv_snapbase = "/n/hernquistfs1/Illustris/SmallBox/GFM/Production/Cosmo/" 
        sb_snapbase = "/n/hernquistfs1/spb/Cosmo/"
        js_snapbase = "/n/home04/jsuresh/runs/Cosmo5_V6/output"
        gam_snapbase = "/n/hernquistfs1/jsuresh/Runs/"

        bmap = brewer2mpl.get_map('Spectral','Diverging',9, reverse=True)
        cmap = bmap.mpl_colormap

        # Cosmo0_V6 (512^3 particles):
        if c0_256:
            for redshift in redshifts:
                self.run_list.append("c0_256")
                self.color_list.append(cmap(0.))
                self.snapdir_list.append(mv_snapbase+"Cosmo0_V6/L25n256/output/")
                if redshift == '4':
                    self.label_list.append("Cold Winds+BHs-256 (z=4)")
                    self.snapnum_list.append(54)
                elif redshift == '3':
                    self.label_list.append("Cold Winds+BHs-256 (z=3)")
                    self.snapnum_list.append(60)
                elif redshift == '2.5':
                    self.label_list.append("Cold Winds+BHs-256 (z=2.5)")
                    self.snapnum_list.append(63)
                elif redshift == '2':
                    # self.label_list.append("Cold Winds+BHs-256 (z=2)")
                    self.label_list.append("Fiducial")
                    self.snapnum_list.append(68)
                elif redshift == '1':
                    raise Exception('z=1 not implemented')
                    self.label_list.append("Cold Winds+BHs-256 (z=1)")
                    self.snapnum_list.append()
                elif redshift == '0.3':
                    self.label_list.append("Cold Winds+BHs-256 (z=0.3)")
                    self.snapnum_list.append(114)
                elif redshift == '0':
                    self.label_list.append("Cold Winds+BHs-256 (z=0)")
                    self.snapnum_list.append(135)
        if c0_512:
            for redshift in redshifts:
                self.run_list.append("c0_512")
                self.color_list.append("black")
                self.snapdir_list.append(mv_snapbase+"Cosmo0_V6/L25n512/output/")
                if redshift == '4':
                    self.label_list.append("Cold Winds+BHs (z=4)")
                    self.snapnum_list.append(54)
                elif redshift == '3':
                    self.label_list.append("Cold Winds+BHs (z=3)")
                    self.snapnum_list.append(60)
                elif redshift == '2':
                    self.label_list.append("Cold Winds+BHs (z=2)")
                    self.snapnum_list.append(68)
                elif redshift == '1':
                    raise Exception('z=1 not implemented')
                    self.label_list.append("Cold Winds+BHs (z=1)")
                    self.snapnum_list.append()
                elif redshift == '0.3':
                    self.label_list.append("Cold Winds+BHs (z=0.3)")
                    self.snapnum_list.append(114)
                elif redshift == '0':
                    self.label_list.append("Cold Winds+BHs (z=0)")
                    self.snapnum_list.append(135)


        # Cosmo2_V6 (256^3 particles):
        if c2_256:
            for redshift in redshifts:
                self.run_list.append("c2_256")
                self.color_list.append("cyan")
                self.snapdir_list.append(mv_snapbase+"Cosmo2_V6/L25n256/output/")
                self.linestyle_list.append("solid")
                if redshift == '4':
                    self.label_list.append("Cold Winds only (z=4)")
                    self.snapnum_list.append(54)
                elif redshift == '3':
                    self.label_list.append("Cold Winds only (z=3)")
                    self.snapnum_list.append(60)
                elif redshift == '2':
                    # self.label_list.append("Cold Winds only (z=2)")
                    self.label_list.append("No AGN")
                    self.snapnum_list.append(68)
                elif redshift == '0.3':
                    self.label_list.append("Cold Winds only (z=0.3)")
                    self.snapnum_list.append(114)
                elif redshift == '0':
                    self.label_list.append("Cold Winds only (z=0)")
                    self.snapnum_list.append(135)
                

        # Cosmo3_V6 (512^3 particles):
        if c3_512:
            for redshift in redshifts:
                self.run_list.append("c3_512")
                self.color_list.append("black")
                self.snapdir_list.append(sb_snapbase+"Cosmo3_V6/L25n512/output/")
                if redshift == '4':
                    self.label_list.append("No Feedback (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("No Feedback (z=3)")
                    self.snapnum_list.append(3)
                elif redshift == '2':
                    self.label_list.append("No Feedback (z=2)")
                    self.snapnum_list.append(5)


        # Cosmo4_V6 (512^3 particles):
        if c4_512:
            for redshift in redshifts:
                self.run_list.append("c4_512")
                self.color_list.append("green")
                self.linestyle_list.append('solid')
                self.snapdir_list.append(sb_snapbase+"Cosmo4_V6/L25n512/output/")
                if redshift == '4':
                    self.label_list.append("Hot Winds Only (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Hot Winds Only (z=4)")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Hot Winds Only (z=2.5)")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    self.label_list.append("Hot Winds Only (z=4)")
                    self.snapnum_list.append(5)
                elif redshift == '0.3':
                    self.label_list.append("Hot Winds Only (z=4)")
                    self.snapnum_list.append(9)

        if c4_check:
            for redshift in redshifts:
                self.run_list.append("c4_check")
                self.color_list.append("red")
                self.linestyle_list.append('solid')
                self.snapdir_list.append(gam_snapbase+"c4_check/output/")
                if redshift == '2.5':
                    self.label_list.append("c4_check")
                    self.snapnum_list.append(4)


        # Cosmo0_V6_fastWinds (256^3 particles):
        if c0_fw_256:
            for redshift in redshifts:
                self.run_list.append("c0_fw_256")
                self.color_list.append("orange")
                self.snapdir_list.append(mv_snapbase+"Cosmo0_V6_fastWinds/L25n256/output/")
                if redshift == '4':
                    self.label_list.append("Fast Winds+BHs (z=4)")
                    self.snapnum_list.append(54)
                elif redshift == '3':
                    self.label_list.append("Fast Winds+BHs (z=3)")
                    self.snapnum_list.append(60)
                elif redshift == '2':
                    self.label_list.append("Fast Winds+BHs (z=2)")
                    self.snapnum_list.append(68)
                elif redshift == '0.3':
                    self.label_list.append("Fast Winds+BHs (z=0.3)")
                    self.snapnum_list.append(114)

        # Cosmo0_V6_strongWinds (256^3 particles):
        if c0_sw_256:
            for redshift in redshifts:
                self.run_list.append("c0_sw_256")
                self.color_list.append("chartreuse")
                self.snapdir_list.append(mv_snapbase+"Cosmo0_V6_strongWinds/L25n256/output/")
                self.linestyle_list.append("solid")
                if redshift == '4':
                    self.label_list.append("Strong Winds+BHs (z=4)")
                    self.snapnum_list.append(54)
                elif redshift == '3':
                    self.label_list.append("Strong Winds+BHs (z=4)")
                    self.snapnum_list.append(60)
                elif redshift == '2':
                    self.label_list.append("Higher Mass-loading")
                    self.snapnum_list.append(68)
                elif redshift == '0.3':
                    self.label_list.append("Strong Winds+BHs (z=4)")
                    self.snapnum_list.append(114)

        # Cosmo5_V6 (256^3 particles):
        if c5_256:
            for redshift in redshifts:
                self.run_list.append("c5_256")
                self.color_list.append("magenta")
                self.snapdir_list.append(js_snapbase)
                if redshift == '4':
                    self.label_list.append("Intermediate-speed Winds Only (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Intermediate-speed Winds Only (z=3)")
                    self.snapnum_list.append(3)
                elif redshift == '2':
                    self.label_list.append("Intermediate-speed Winds Only (z=2)")
                    self.snapnum_list.append(5)
            
        # gamma boxes (256^3 particles): 
        if g0_BH:
            for redshift in redshifts:
                self.run_list.append("c0_256")
                #self.color_list.append(cmap(0.))
                self.color_list.append('blue')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(mv_snapbase+"Cosmo0_V6/L25n256/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0+BHs (z=4)")
                    self.snapnum_list.append(54)
                elif redshift == '3':
                    self.label_list.append("Gamma=0+BHs (z=3)")
                    self.snapnum_list.append(60)
                elif redshift == '2.5':
                    self.label_list.append("Gamma=0+BHs (z=2.5)")
                    self.snapnum_list.append(63)
                elif redshift == '2':
                    # self.label_list.append("Gamma=0+BHs (z=2)")
                    self.label_list.append("Fiducial")
                    self.snapnum_list.append(68)
                elif redshift == '1':
                    raise Exception('z=1 not implemented')
                    self.label_list.append("Gamma=0+BHs (z=1)")
                    self.snapnum_list.append()
                elif redshift == '0.3':
                    self.label_list.append("Gamma=0+BHs (z=0.3)")
                    self.snapnum_list.append(114)
                elif redshift == '0':
                    self.label_list.append("Gamma=0+BHs (z=0)")
                    self.snapnum_list.append(135)

        if g10_BH:
            for redshift in redshifts:
                self.run_list.append("gam_10_BH")
                self.color_list.append(cmap(0.2))
                self.linestyle_list.append('solid')
                self.snapdir_list.append(gam_snapbase+"gam_10_BH/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0.1 + BHs (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Gamma=0.1 + BHs (z=3)")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Gamma=0.1 + BHs (z=2.5)")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    self.label_list.append("Gamma=0.1 + BHs (z=2)")
                    self.snapnum_list.append(5)
                elif redshift == '1':
                    self.label_list.append("Gamma=0.1 + BHs (z=1)")
                    self.snapnum_list.append(7)
                elif redshift == '0.3':
                    self.label_list.append("Gamma=0.1 + BHs (z=0.3)")
                    self.snapnum_list.append(9)

        if g20_BH:
            for redshift in redshifts:
                self.run_list.append("gam_20_BH")
                self.linestyle_list.append('solid')
                self.color_list.append(cmap(0.4))
                self.snapdir_list.append(gam_snapbase+"gam_20_BH/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0.2 + BHs (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Gamma=0.2 + BHs (z=3)")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Gamma=0.2 + BHs (z=2.5)")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    self.label_list.append("Gamma=0.2 + BHs (z=2)")
                    self.snapnum_list.append(5)
                elif redshift == '1':
                    self.label_list.append("Gamma=0.2 + BHs (z=1)")
                    self.snapnum_list.append(7)
                elif redshift == '0.3':
                    self.label_list.append("Gamma=0.2 + BHs (z=0.3)")
                    self.snapnum_list.append(9)

        if g25_BH:
            for redshift in redshifts:
                self.run_list.append("gam_25_BH")
                self.linestyle_list.append('solid')
                self.color_list.append("blue")
                self.snapdir_list.append(gam_snapbase+"gam_25_BH/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0.25 + BHs (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Gamma=0.25 + BHs (z=1)")
                    self.snapnum_list.append(3)
                elif redshift == '2':
                    self.label_list.append("Gamma=0.25 + BHs (z=0.3)")
                    self.snapnum_list.append(5)

        if g30_BH:
            for redshift in redshifts:
                self.run_list.append("gam_30_BH")
                self.linestyle_list.append('solid')
                self.color_list.append(cmap(0.6))
                self.snapdir_list.append(gam_snapbase+"gam_30_BH/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0.3 + BHs (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Gamma=0.3 + BHs (z=3)")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Gamma=0.3 + BHs (z=2.5)")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    self.label_list.append("Gamma=0.3 + BHs (z=2)")
                    self.snapnum_list.append(5)
                elif redshift == '1':
                    self.label_list.append("Gamma=0.3 + BHs (z=1)")
                    self.snapnum_list.append(7)
                elif redshift == '0.3':
                    self.label_list.append("Gamma=0.3 + BHs (z=0.3)")
                    self.snapnum_list.append(9)

        if g40_BH:
            for redshift in redshifts:
                self.run_list.append("gam_40_BH")
                self.linestyle_list.append('solid')
                self.color_list.append(cmap(0.8))
                self.snapdir_list.append(gam_snapbase+"gam_40_BH/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0.4 + BHs (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Gamma=0.4 + BHs (z=3)")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Gamma=0.4 + BHs (z=2.5)")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    self.label_list.append("Gamma=0.4 + BHs (z=2)")
                    self.snapnum_list.append(5)
                elif redshift == '1':
                    self.label_list.append("Gamma=0.4 + BHs (z=1)")
                    self.snapnum_list.append(7)
                elif redshift == '0.3':
                    self.label_list.append("Gamma=0.4 + BHs (z=0.3)")
                    self.snapnum_list.append(9)

        if g50_BH:
            for redshift in redshifts:
                self.run_list.append("gam_50_BH")
                self.linestyle_list.append('solid')
                # self.color_list.append(cmap(1.0))
                self.color_list.append('saddlebrown')
                self.snapdir_list.append(gam_snapbase+"gam_50_BH/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0.5 + BHs (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Gamma=0.5 + BHs (z=3)")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Gamma=0.5 + BHs (z=2.5)")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    # self.label_list.append("Gamma=0.5 + BHs (z=2)")
                    self.label_list.append("Fixed-E Hot Winds")
                    self.snapnum_list.append(5)
                elif redshift == '1':
                    self.label_list.append("Gamma=0.5 + BHs (z=1)")
                    self.snapnum_list.append(7)
                elif redshift == '0.3':
                    self.label_list.append("Gamma=0.5 + BHs (z=0.3)")
                    self.snapnum_list.append(9)

        if g75_BH:
            for redshift in redshifts:
                self.run_list.append("gam_75_BH")
                self.linestyle_list.append('solid')
                self.color_list.append("orange")
                self.snapdir_list.append(gam_snapbase+"gam_75_BH/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0.75 + BHs (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Gamma=0.75 + BHs (z=3)")
                    self.snapnum_list.append(3)
                elif redshift == '2':
                    self.label_list.append("Gamma=0.75 + BHs (z=2)")
                    self.snapnum_list.append(5)

        if g95_BH:
            for redshift in redshifts:
                self.run_list.append("gam_95_BH")
                self.linestyle_list.append('solid')
                self.color_list.append("red")
                self.snapdir_list.append(gam_snapbase+"gam_95_BH/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0.95 + BHs (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Gamma=0.95 + BHs (z=3)")
                    self.snapnum_list.append(3)
                elif redshift == '2':
                    self.label_list.append("Gamma=0.95 + BHs (z=2)")
                    self.snapnum_list.append(5)

        # No black hole runs:
        if g10_noBH:
            for redshift in redshifts:
                self.run_list.append("gam_10_noBH")
                self.linestyle_list.append('dashed')
                self.color_list.append(cmap(0.2))
                self.snapdir_list.append(gam_snapbase+"gam_10_noBH/output/")
                if redshift == '2.5':
                    self.label_list.append(None)
                    self.snapnum_list.append(4)
                if redshift == '2':
                    self.label_list.append(None)
                    self.snapnum_list.append(5)
                elif redshift == '0.3':
                    self.label_list.append("Gamma=0.1-noBH (z=0.3)")
                    self.snapnum_list.append(9)

        if g20_noBH:
            for redshift in redshifts:
                self.run_list.append("gam_20_noBH")
                self.linestyle_list.append('dashed')
                self.color_list.append(cmap(0.4))
                self.snapdir_list.append(gam_snapbase+"gam_20_noBH/output/")
                if redshift == '2':
                    self.label_list.append(None)
                    self.snapnum_list.append(5)

        if g30_noBH:
            for redshift in redshifts:
                self.run_list.append("gam_30_noBH")
                self.linestyle_list.append('dashed')
                self.color_list.append(cmap(0.6))
                self.snapdir_list.append(gam_snapbase+"gam_30_noBH/output/")
                if redshift == '2':
                    self.label_list.append(None)
                    self.snapnum_list.append(5)

        if g40_noBH:
            for redshift in redshifts:
                self.run_list.append("gam_40_noBH")
                self.linestyle_list.append('dashed')
                self.color_list.append(cmap(0.8))
                self.snapdir_list.append(gam_snapbase+"gam_40_noBH/output/")
                if redshift == '2':
                    self.label_list.append(None)
                    self.snapnum_list.append(5)


        if g50_noBH:
            for redshift in redshifts:
                self.run_list.append("gam_50_noBH")
                self.linestyle_list.append('dashed')
                self.color_list.append(cmap(1.0))
                self.snapdir_list.append(gam_snapbase+"gam_50_noBH/output/")
                if redshift == '2.5':
                    self.label_list.append(None)
                    self.snapnum_list.append(4)
                if redshift == '2':
                    self.label_list.append(None)
                    self.snapnum_list.append(5)
                elif redshift == '0.3':
                    self.label_list.append("Gamma=0.5-noBH (z=0.3)")
                    self.snapnum_list.append(9)

        # No thermal component runs:
        if g10_nothermal:
            for redshift in redshifts:
                self.run_list.append("gam_10_BH_nothermal")
                self.linestyle_list.append('dashed')
                self.color_list.append(cmap(0.2))
                self.snapdir_list.append(gam_snapbase+"gam_10_BH_nothermal/output/")
                if redshift == '2.5':
                    self.label_list.append(None)
                    self.snapnum_list.append(4)
                if redshift == '2':
                    self.label_list.append(None)
                    self.snapnum_list.append(5)
                elif redshift == '0.3':
                    self.label_list.append("Gamma=0.1+BHs-thermal (z=0.3)")
                    self.snapnum_list.append(9)

        if g20_nothermal:
            for redshift in redshifts:
                self.run_list.append("gam_20_BH_nothermal")
                self.linestyle_list.append('dashed')
                self.color_list.append(cmap(0.4))
                self.snapdir_list.append(gam_snapbase+"gam_20_BH_nothermal/output/")
                if redshift == '2':
                    self.label_list.append(None)
                    self.snapnum_list.append(5)

        if g30_nothermal:
            for redshift in redshifts:
                self.run_list.append("gam_30_BH_nothermal")
                self.linestyle_list.append('dashed')
                self.color_list.append(cmap(0.6))
                self.snapdir_list.append(gam_snapbase+"gam_30_BH_nothermal/output/")
                if redshift == '2':
                    self.label_list.append(None)
                    self.snapnum_list.append(5)

        if g40_nothermal:
            for redshift in redshifts:
                self.run_list.append("gam_40_BH_nothermal")
                self.linestyle_list.append('dashed')
                self.color_list.append(cmap(0.8))
                self.snapdir_list.append(gam_snapbase+"gam_40_BH_nothermal/output/")
                if redshift == '2':
                    self.label_list.append(None)
                    self.snapnum_list.append(5)

        if g50_nothermal:
            for redshift in redshifts:
                self.run_list.append("gam_50_BH_nothermal")
                self.linestyle_list.append('dashed')
                self.color_list.append(cmap(1.0))
                self.snapdir_list.append(gam_snapbase+"gam_50_BH_nothermal/output/")
                if redshift == '2.5':
                    self.label_list.append(None)
                    self.snapnum_list.append(4)
                if redshift == '2':
                    self.label_list.append(None)
                    self.snapnum_list.append(5)
                elif redshift == '0.3':
                    self.label_list.append("Gamma=0.5+BHs-thermal (z=0.3)")
                    self.snapnum_list.append(9)


        # Class I wind runs (fixv):
        if g25_fixv:
            for redshift in redshifts:
                self.run_list.append("gam_25_fixv")
                self.color_list.append("orange")
                self.linestyle_list.append("solid")
                self.snapdir_list.append(gam_snapbase+"gam_25_fixv/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0.25-fixv (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Gamma=0.25-fixv (z=3)")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Gamma=0.25-fixv (z=2.5)")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    self.label_list.append("Gamma=0.25-fixv (z=2)")
                    self.snapnum_list.append(5)

        if g50_fixv:
            for redshift in redshifts:
                self.run_list.append("gam_50_fixv")
                self.color_list.append("magenta")
                self.linestyle_list.append("solid")
                self.snapdir_list.append(gam_snapbase+"gam_50_fixv/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0.5-fixv (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Gamma=0.5-fixv (z=3)")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Gamma=0.5-fixv (z=2.5)")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    # self.label_list.append("Gamma=0.5-fixv (z=2)")
                    self.label_list.append("Fixed-v Hot Winds")
                    self.snapnum_list.append(5)
                elif redshift == '0.3':
                    self.label_list.append("Gamma=0.5-fixv (z=0.3)")
                    self.snapnum_list.append(9)

        if g25_fixv_nothermal:
            for redshift in redshifts:
                self.run_list.append("gam_25_fixv_nothermal")
                self.color_list.append("orange")
                self.linestyle_list.append("dashed")
                self.snapdir_list.append(gam_snapbase+"gam_25_fixv_nothermal/output/")
                if redshift == '2':
                    self.label_list.append(None)
                    self.snapnum_list.append(5)

        if g50_fixv_nothermal:
            for redshift in redshifts:
                self.run_list.append("gam_50_fixv_nothermal")
                self.color_list.append("gold")
                self.linestyle_list.append("solid")
                self.snapdir_list.append(gam_snapbase+"gam_50_fixv_nothermal/output/")
                if redshift == '2':
                    self.label_list.append("Faster Winds")
                    self.snapnum_list.append(5)


        # Class III wind runs (fixv):
        if g25_fixv_fixeta:
            for redshift in redshifts:
                self.run_list.append("gam_25_fixv_fixeta")
                self.color_list.append("orange")
                self.linestyle_list.append("solid")
                self.snapdir_list.append(gam_snapbase+"gam_25_fixv_fixeta/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0.25-fixv_fixeta (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Gamma=0.25-fixv_fixeta (z=3)")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Gamma=0.25-fixv_fixeta (z=2.5)")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    self.label_list.append("Gamma=0.25-fixv_fixeta (z=2)")
                    self.snapnum_list.append(5)

        if g50_fixv_fixeta:
            for redshift in redshifts:
                self.run_list.append("gam_50_fixv_fixeta")
                self.color_list.append("red")
                self.linestyle_list.append("solid")
                self.snapdir_list.append(gam_snapbase+"gam_50_fixv_fixeta/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0.5-fixv_fixeta (z=4)")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Gamma=0.5-fixv_fixeta (z=3)")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Gamma=0.5-fixv_fixeta (z=2.5)")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    self.label_list.append("Gamma=0.5-fixv_fixeta (z=2)")
                    self.snapnum_list.append(5)

        # Wind Enrichment:
        if c0_nometalwinds:
            for redshift in redshifts:
                self.run_list.append("c0_nometalwinds")
                self.color_list.append("gray")
                self.linestyle_list.append("solid")
                self.snapdir_list.append(gam_snapbase+"c0_nometalwinds/output/")
                if redshift == '4':
                    self.label_list.append("Pristine Winds")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Pristine Winds")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Pristine Winds")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    self.label_list.append("Pristine Winds")
                    self.snapnum_list.append(5)

        if c0_fullmetalwinds:
            for redshift in redshifts:
                self.run_list.append("c0_fullmetalwinds")
                self.color_list.append("black")
                self.linestyle_list.append("solid")
                self.snapdir_list.append(gam_snapbase+"c0_fullmetalwinds/output/")
                if redshift == '4':
                    self.label_list.append("Fully Enriched Winds")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("Fully Enriched Winds")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Fully Enriched Winds")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    self.label_list.append("Fully Enriched Winds")
                    self.snapnum_list.append(5)

        self.ndat = len(self.run_list)

if __name__ == '__main__':
    full_analysis()

