import os
import time
import numpy as np
import h5py

import matplotlib
matplotlib.use('PDF')
from matplotlib import rc
# rc('text', usetex=True)
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
        # self.CGMsnap_base = '/n/home04/jsuresh/CGM_new/data/CGM_snaps/'
        # self.grid_base = '/n/home04/jsuresh/CGM_new/data/grids/'
        # self.fig_base = '/n/home04/jsuresh/CGM_new/data/figs/'
        # self.npz_base = '/n/home04/jsuresh/CGM_new/data/npz/'
        self.CGMsnap_base = '/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/CGM_new/data/CGM_snaps/'
        self.grid_base = '/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/CGM_new/data/grids/'
        self.fig_base = '/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/CGM_new/data/figs/'
        self.npz_base = '/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/CGM_new/data/npz/'

        # Select which runs to be analyzed concurrently:
        self.run_list = []
        self.color_list = [] 
        self.label_list = [] 
        self.linestyle_list = []
        self.snapnum_list = [] 
        self.snapdir_list = []


        # Feedback paper runs:
        # redshifts = ['2.5']
        # self.dat_prep(redshifts,g0_BH=1)
        # self.dat_prep(redshifts,c2_256=1)
        # self.dat_prep(redshifts,c0_sw_256=1)
        # self.dat_prep(redshifts,g50_fixv_nothermal=1)
        # self.dat_prep(redshifts,g50_BH=1)
        # self.dat_prep(redshifts,g50_fixv=1)
        # self.dat_prep(redshifts,c0_nometalwinds=1)
        # self.dat_prep(redshifts,c0_fullmetalwinds=1)
        # self.dat_prep(redshifts,c0_sw_256=1)

        # Check runs:
        redshifts = ['2.5']
        # self.dat_prep(redshifts,c2_256=1)
        # self.dat_prep(redshifts,c0_128=1)
        # self.dat_prep(redshifts,g0_BH=1)
        # self.dat_prep(redshifts,c0_512=1)
        # self.dat_prep(redshifts,g0_BH=1)
        # self.dat_prep(redshifts,c0_check=1)
        # self.dat_prep(redshifts,c0_nometalwinds=1)
        # self.dat_prep(redshifts,c0_fullmetalwinds=1)
        # self.dat_prep(redshifts,c0_dumpfactor95=1)
        # self.dat_prep(redshifts,no_metal_cooling=1)
        # self.dat_prep(redshifts,low_recouple_dens=1)
        # self.dat_prep(redshifts,high_recouple_dens=1)
        self.dat_prep(redshifts,multi_speed=1)
        # self.dat_prep(redshifts,g50_fixv_nothermal=1)



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
        # self.phase_budget('m',mmin=10.**11,savename='m_budget_res',minimal=True)
        # self.phase_budget('z',mmin=10.**11,savename='z_budget_res_recouple',minimal=True)
        # self.CGM_and_gal_metallicity(savename="massmet_newleg")
        # self.radial_profile('T',savename='T_res')
        # self.radial_profile('z',savename='z_res')
        # self.columndens_vs_R("O6",150.,minmass=10**11.8,maxmass=10**12.2,plot_mode=1,savename="O6_z3_galaxybinned")

        ##############################
        # 2D data analysis functions #
        ##############################
        # self.plot_grids("H1")
        # self.grid_sightlines("H1",200.,coldens_min=15.5,minmass=10**11.8,maxmass=10**12.2,coverfrac_within_R=True,rudie_155=True,savename="rudie_155_z2_recouple",show_Fumagalli=True)
        # self.grid_sightlines("H1",200.,coldens_min=17.2,minmass=10**11.8,maxmass=10**12.2,coverfrac_within_R=True,rudie_172=True,savename="rudie_172_z2_recouple",show_Fumagalli=True)
        # self.grid_sightlines("H1",200.,coldens_min=19.,minmass=10**11.8,maxmass=10**12.2,coverfrac_within_R=True,rudie_19=True,savename="rudie_19_z2_recouple",show_Fumagalli=True)
        # self.grid_sightlines("H1",200.,coldens_min=20.3,minmass=10**11.8,maxmass=10**12.2,coverfrac_within_R=True,rudie_203=True,savename="rudie_203_z2_recouple",show_Fumagalli=True)
        # print "2"
        # self.grid_sightlines("H1",200.,coldens_min=20.3,minmass=10**11.8,maxmass=10**12.2,coverfrac_within_R=True,rudie_DLA=True,savename="DLA_res")
        # print "3"
        # self.grid_sightlines("H1",200.,minmass=10**11.8,maxmass=10**12.2,coldens_vs_R=True,savename="coldens_newfan2")
        # self.grid_sightlines("O6",150.,minmass=10**11.8,maxmass=10**12.2,coldens_vs_R=True,savename="O6_z25_coldens")

        # self.grid_sightlines("H1",300.,minmass=10**11.8,maxmass=10**12.2,coldens_vs_R=True,savename="H1_z25_coldens_300")
        # self.grid_sightlines("Mg2",300.,minmass=10**11.8,maxmass=10**12.2,coldens_vs_R=True,savename="Mg2_z25_coldens_300")
        # self.grid_sightlines("C3",300.,minmass=10**11.8,maxmass=10**12.2,coldens_vs_R=True,savename="C3_z25_coldens_300")
        # self.grid_sightlines("C4",300.,minmass=10**11.8,maxmass=10**12.2,coldens_vs_R=True,turner_c4=True,savename="C4_z2_multispeed")
        # self.grid_sightlines("Si3",300.,minmass=10**11.8,maxmass=10**12.2,coldens_vs_R=True,savename="Si3_z25_coldens_300")
        # self.grid_sightlines("Si4",300.,minmass=10**11.8,maxmass=10**12.2,coldens_vs_R=True,savename="Si4_z25_coldens_300")
        # self.grid_sightlines("O6",300.,minmass=10**11.8,maxmass=10**12.2,coldens_vs_R=True,savename="O6_z2_coldens_300_res")
        self.grid_sightlines("C-1",300.,minmass=10**11.8,maxmass=10**12.2,coldens_vs_R=True,turner_c4=True,savename="C_z2_multispeed")

        # self.plot_grids("H1",HI_custom_map=True)
        # self.plot_grids("C4",vmax=15.,vmin=5.)
        # self.plot_grids("C-1") #,vmax=17.,vmin=5.
        # self.plot_grids("Si4")
        # self.plot_grids("O6")
        # self.coverfrac_bootstrap("H1",200.,coldens_min=17.2,minmass=10**11.8,maxmass=10**12.2)
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
                p1.semilogy(np.log10(mbins_med),med[:,0]/((0.0456/0.27)*mbins_med),color=self.color_list[i],linestyle='dashed')
                p1.semilogy(np.log10(mbins_med),med[:,1]/((0.0456/0.27)*mbins_med),color=self.color_list[i],label=self.label_list[i])
                # legend1 = p1.legend([p_gal,p_CGM], ["Stars+ISM","CGM"], prop={'size':6})
                # plt.gca().add_artist(legend1)
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
        p1.set_xlabel(r"$\log_{10}\left[\frac{M_{\rm Halo}}{M_\odot}\right]$")
        # p1.set_ylabel(r"$\log_{10}\left[ \frac{M_{\rm baryons}}{M_{\rm Halo} \times \Omega_b} \right]$")
        p1.set_ylabel(r"$\frac{M_{\rm baryons}}{M_{\rm Halo} \times f_b}$")
        # lg = p1.legend(loc=2,prop={'size':5.5},ncol=2,columnspacing=0.3,borderpad=0.2)
        lg = p1.legend(loc=2,prop={'size':6.1},ncol=2,columnspacing=0.8,handletextpad=0.1,handlelength=1.6)
        lg.draw_frame(False)

        font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 14,
        }
        p1.text(11.05,2.4,"CGM",fontdict=font)
        p1.text(11.05,0.3,"ISM+Stars",fontdict=font)
        p2.text(12.,0.02,"CGM",fontdict=font)
        p2.text(11.05,0.6,"ISM+Stars",fontdict=font)
        # plt.gca().add_artist(legend1)

        p2.xaxis.set_ticks(np.arange(11.0,12.6,0.5))
        # p2.set_ylabel(r"$\log_{10}\left[ \frac{Z_{\rm baryons}}{Z_\odot} \right]$")
        p2.set_ylabel(r"$\frac{Z_{\rm baryons}}{Z_\odot}$")
        p2.set_xlabel(r"$\log_{10}\left[\frac{M_{\rm Halo}}{M_\odot}\right]$")
        #p2.legend(loc=2,prop={'size':5})

        plt.subplots_adjust(wspace=0.4)

        if savename != None: 
            plt.savefig(self.fig_base+savename+".pdf", bbox_inches='tight')
        else:
            plt.savefig(self.fig_base+"mass_budget.pdf", bbox_inches='tight')




    def phase_budget(self,data_type,mmin=1e11,mmax=1e18,savename=None,log_plot=False,minimal=True):
        plt.close('all')    
        if not minimal:
            plt.figure(figsize=(3.54331,3.14*2))
            p1 = plt.subplot(4,1,1)
            p2 = plt.subplot(4,1,2)
            p3 = plt.subplot(4,1,3)
            p4 = plt.subplot(4,1,4)
        else:
            plt.figure(figsize=(3.54331,3.14*2))
            p1 = plt.subplot(2,1,1)
            p2 = plt.subplot(2,1,2)

        bmap = brewer2mpl.get_map('Spectral','Diverging',9, reverse=True)
        cmap = bmap.mpl_colormap

        rho_cut = 10**(-2.8)
        rho_SFR = 0.13
        T_cut = 10.**5.

        # Bin by galaxy stellar mass:
        #mbins_min = 10.**(9.+np.array([0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0,1.2,1.4,1.6,1.8]))
        #mbins_max = 10.**(9.+np.array([.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0,1.2,1.4,1.6,1.8,2.0]))
        mbins_min = 10.**(9.+np.array([.5,.6,.7,.8,.9,1.0,1.2,1.4,1.6,1.8]))
        mbins_max = 10.**(9.+np.array([.6,.7,.8,.9,1.0,1.2,1.4,1.6,1.8,2.0]))
        mbins_med = 10.**((np.log10(mbins_min) + np.log10(mbins_max))/2.)
        n_mbins = np.size(mbins_min)

        for i in np.arange(self.ndat):
            npz_fname = "phasebudget_z_{}_s{}_9mingalmass.npz".format(self.run_list[i],self.snapnum_list[i])
            if os.path.isfile(self.npz_base+npz_fname):
                print "Loaded from npz!: {}".format(self.npz_base+npz_fname)
                f = np.load(self.npz_base+npz_fname)
                gal_mass = f['gal_mass']
                vw_budget_gal = f['vw_budget_gal']
                mw_budget_gal = f['mw_budget_gal']
                zw_budget_gal = f['zw_budget_gal']
                f.close()
            else:
                [grp_ids,grp_mass] = self.find_desired_grpids(i,minmass=mmin,maxmass=mmax,return_grpmass=True)
                gal_dict = self.get_gal_props(i,grp_ids)
                gal_mass = gal_dict['gal_mass']
                in_mbin = gal_mass > 10.**9.5
                gal_mass = gal_mass[in_mbin]
                selected_grp_ids = grp_ids[in_mbin]
                # selected_grp_ids = grp_ids
                # gal_mass = grp_mass/100. # clearly wrong - just for plot sake

                # no radial cut implemented yet!
                # m = np.zeros(0)
                # z = np.zeros(0)
                # T = np.zeros(0)
                # rho = np.zeros(0)
                # vol = np.zeros(0)
                

                mw_budget_gal = np.zeros([np.size(selected_grp_ids),4])
                zw_budget_gal = np.zeros([np.size(selected_grp_ids),4])
                vw_budget_gal = np.zeros([np.size(selected_grp_ids),4])
                for j in np.arange(np.size(selected_grp_ids)):
                    print "{} of {}".format(j,np.size(selected_grp_ids))
                    grp_id = selected_grp_ids[j]
                    grp_data = self.load_CGMsnap_data(i,grp_id,data_entry='all')
                    redshift = grp_data['redshift']
                    #pos = AU.PhysicalPosition(np.array(grp_data['Coordinates']),grp_data['redshift'])
                    #grp_pos = AU.PhysicalPosition(np.array(grp_data['grp_pos']),grp_data['redshift'])
                    #r_los = np.sqrt((pos[:,0]-grp_pos[0])**2+(pos[:,1]-grp_pos[1])**2)
                    m = AU.PhysicalMass(np.array(grp_data['Masses']))
                    z = np.array(grp_data['GFM_Metallicity'])
                    u = np.array(grp_data['InternalEnergy'])
                    nelec = np.array(grp_data["ElectronAbundance"])
                    T = AU.GetTemp(u, nelec, gamma = 5.0/3.0)
                    rho = AU.PhysicalDensity(np.array(grp_data['Density']),redshift)
                    vol = AU.PhysicalVolume(grp_data['Volume'],redshift)

                    rho /= AU.ProtonMass

                    cde = np.logical_and(np.logical_and(T < T_cut,rho >= rho_cut),rho < rho_SFR)
                    cdi = np.logical_and(T < T_cut,rho < rho_cut)
                    hde = np.logical_and(np.logical_and(T > T_cut,rho >= rho_cut),rho < rho_SFR)
                    hdi = np.logical_and(T > T_cut,rho < rho_cut)
                    CGM = rho < rho_SFR
                    m_tot = np.sum(m[CGM])
                    z_tot = np.sum(m[CGM]*z[CGM])
                    v_tot = np.sum(vol[CGM])
                    mw_budget = np.array([np.sum(m[cde]),np.sum(m[cdi]),np.sum(m[hde]),np.sum(m[hdi])])/m_tot
                    zw_budget = np.array([np.sum(m[cde]*z[cde]),np.sum(m[cdi]*z[cdi]),np.sum(m[hde]*z[hde]),np.sum(m[hdi]*z[hdi])])/z_tot
                    # print "np.sum(cde) ",np.sum(cde)
                    # print "np.min(T) ",np.min(T)
                    # print "np.max(rho) ",np.min(T)
                    vw_budget = np.array([np.sum(vol[cde]),np.sum(vol[cdi]),np.sum(vol[hde]),np.sum(vol[hdi])])/v_tot
                    mw_budget_gal[j,:] = mw_budget
                    zw_budget_gal[j,:] = zw_budget
                    vw_budget_gal[j,:] = vw_budget


                # np.savez(self.npz_base+npz_fname,gal_mass=gal_mass,mw_budget_gal=mw_budget_gal,zw_budget_gal=zw_budget_gal,vw_budget_gal=vw_budget_gal)


            # if data_type == 'v':
            #     [Q1,med,Q3] = AU._calc_percentiles(vw_budget_gal)
            # elif data_type == 'm':
            #     [Q1,med,Q3] = AU._calc_percentiles(mw_budget_gal)
            # elif data_type == 'z':
            #     [Q1,med,Q3] = AU._calc_percentiles(zw_budget_gal)

            if data_type == 'v':
                u = vw_budget_gal
                p1.set_title(r'CGM Volume Budget')
            elif data_type == 'm':
                u = mw_budget_gal
                # p1.set_title(r'CGM Mass Budget')
                p1.set_ylabel(r'Mass Fraction')
                p2.set_ylabel(r'Mass Fraction')
            elif data_type == 'z':
                u = zw_budget_gal
                # p1.set_title(r'CGM Metal Budget')
                p1.set_ylabel(r'Metal Fraction')
                p2.set_ylabel(r'Metal Fraction')

            # print "vw_budget_gal ",vw_budget_gal
            # print "mw_budget_gal ",mw_budget_gal
            # print "zw_budget_gal ",zw_budget_gal

            [Q1_cde,med_cde,Q3_cde] = AU._calc_percentiles_v2(gal_mass,mbins_min,mbins_max,u[:,0])
            [Q1_cdi,med_cdi,Q3_cdi] = AU._calc_percentiles_v2(gal_mass,mbins_min,mbins_max,u[:,1])
            [Q1_hde,med_hde,Q3_hde] = AU._calc_percentiles_v2(gal_mass,mbins_min,mbins_max,u[:,2])
            [Q1_hdi,med_hdi,Q3_hdi] = AU._calc_percentiles_v2(gal_mass,mbins_min,mbins_max,u[:,3])

            if log_plot:
                p1.loglog(mbins_med,med_cde,label=self.label_list[i],color=self.color_list[i],linestyle=self.linestyle_list[i])
                p1.set_ylim([0.05,1.])
                p1.set_xlim([10.**9.5,10.**11])
                p2.loglog(mbins_med,med_cdi,label=self.label_list[i],color=self.color_list[i],linestyle=self.linestyle_list[i])
                p2.set_ylim([0.05,1.])
                p2.set_xlim([10.**9.5,10.**11])
                p3.loglog(mbins_med,med_hde,label=self.label_list[i],color=self.color_list[i],linestyle=self.linestyle_list[i])
                p3.set_ylim([0.05,1.])
                p3.set_xlim([10.**9.5,10.**11])
                p4.loglog(mbins_med,med_hdi,label=self.label_list[i],color=self.color_list[i],linestyle=self.linestyle_list[i])
                p4.set_ylim([0.05,1.])
                p4.set_xlim([10.**9.5,10.**11])
                p4.set_xlabel('Galaxy Mass')
            else:
                if not minimal: 
                    p1.semilogx(mbins_med,med_cde,label=self.label_list[i],color=self.color_list[i],linestyle=self.linestyle_list[i])
                    p1.set_ylim([0.,1.])
                    p1.set_xlim([10.**9.5,10.**11])
                    p1.set_ylabel(r'Cold-Dense')
                    p2.semilogx(mbins_med,med_cdi,label=self.label_list[i],color=self.color_list[i],linestyle=self.linestyle_list[i])
                    p2.set_ylim([0.,1.])
                    p2.set_xlim([10.**9.5,10.**11])
                    p2.set_ylabel(r'Cold-Diffuse')
                    p3.semilogx(mbins_med,med_hde,label=self.label_list[i],color=self.color_list[i],linestyle=self.linestyle_list[i])
                    p3.set_ylim([0.,1.])
                    p3.set_xlim([10.**9.5,10.**11])
                    p3.set_ylabel(r'Warm-Dense')
                    p4.semilogx(mbins_med,med_hdi,label=self.label_list[i],color=self.color_list[i],linestyle=self.linestyle_list[i])
                    p4.set_ylim([0.,1.])
                    p4.set_xlim([10.**9.5,10.**11])
                    p4.set_ylabel(r'Warm-Diffuse')
                    p4.set_xlabel(r'Galaxy Mass')
                else:
                    font_c = {'family' : 'serif',
                    'color'  : 'darkblue',
                    'weight' : 'normal',
                    'size'   : 14,
                    }

                    font_w = {'family' : 'serif',
                    'color'  : 'darkred',
                    'weight' : 'normal',
                    'size'   : 14,
                    }

                    p1.semilogx(mbins_med,med_cde,label=self.label_list[i],color=self.color_list[i],linestyle=self.linestyle_list[i])
                    p1.set_xlim([10.**9.5,10.**11])
                    p1.set_ylim([0.,1.])
                    # p1.set_ylabel(r'Cool-Dense')
                    if data_type == 'm':
                        p1.text(10.**9.6,0.87,"Cool-Dense",fontdict=font_c)
                    elif data_type == 'z':
                        p1.text(10.**10.9,0.87,"Cool-Dense",fontdict=font_c,ha='right')
                    p1.set_xlabel(r"$\log_{10}\left[\frac{M_{\rm *,gal}}{M_\odot}\right]$")
                    p2.semilogx(mbins_med,med_hdi,label=self.label_list[i],color=self.color_list[i],linestyle=self.linestyle_list[i])
                    p2.set_ylim([0.,1.])
                    p2.set_xlim([10.**9.5,10.**11])
                    # p2.set_ylabel(r'Warm-Diffuse')
                    if data_type == 'm':
                        p2.text(10.**9.6,0.87,"Warm-Diffuse",fontdict=font_w)
                    elif data_type == 'z':
                        p2.text(10.**10.9,0.87,"Warm-Diffuse",fontdict=font_w,ha='right')
                    p2.set_xlabel(r"$\log_{10}\left[\frac{M_{\rm *,gal}}{M_\odot}\right]$")


        if not minimal:
            p2.legend(prop={'size':5})
            plt.subplots_adjust(left=0.2,hspace=0.7)
        else:
            # p2.legend(bbox_to_anchor=(0., 1.05, 1., .107), loc=8, ncol=2, borderaxespad=0.,prop={'size':7})
            # plt.subplots_adjust(left=0.2,hspace=0.6)
            p2.legend(bbox_to_anchor=(0., 2.5, 1., 0.), loc=8, ncol=2, borderaxespad=0.,prop={'size':7})
            plt.subplots_adjust(left=0.2,top=0.8,hspace=0.4)

        if savename != None:
            plt.savefig(self.fig_base+savename+".pdf")
        else:
            plt.savefig(self.fig_base+"{}_budget.pdf".format(data_type))







    def radial_profile(self,data_type,mmin=10**11.8,mmax=10**12.2,savename=None):
        plt.close('all')
        plt.figure(figsize=(3.54331,3.14))

        # Set up bins for radius:
        n_Rbins = 50
        [Rbins_min,Rbins_max] = AU._bin_setup(0.,3.,n_Rbins)
        Rbins_med = (Rbins_min+Rbins_max)/2.

        for i in np.arange(self.ndat):
            npz_fname = "radialprofile_{}_s{}_1012halo.npz".format(self.run_list[i],self.snapnum_list[i])
            if os.path.isfile(self.npz_base+npz_fname):
                print "Loaded from npz!: {}".format(self.npz_base+npz_fname)
                f = np.load(self.npz_base+npz_fname)
                z_profile_gal = f['z_profile_gal']
                T_profile_gal = f['T_profile_gal']
                rho_profile_gal = f['rho_profile_gal']
                f.close()
            else:
                grp_ids = self.find_desired_grpids(i,minmass=mmin,maxmass=mmax)
                # gal_dict = self.get_gal_props(i,grp_ids)
                # gal_mass = gal_dict['gal_mass']
                n_grps = np.size(grp_ids)

                z_profile_gal = np.zeros([n_grps,n_Rbins])
                T_profile_gal = np.zeros([n_grps,n_Rbins])
                rho_profile_gal = np.zeros([n_grps,n_Rbins])
                for j in np.arange(n_grps):
                    print "{} of {}".format(j,n_grps)
                    grp_id = grp_ids[j]
                    grp_data = self.load_CGMsnap_data(i,grp_id,data_entry='all')
                    redshift = grp_data['redshift']
                    r = self._calc_radii(grp_data)
                    m = AU.PhysicalMass(np.array(grp_data['Masses']))
                    z = np.array(grp_data['GFM_Metallicity'])
                    u = np.array(grp_data['InternalEnergy'])
                    nelec = np.array(grp_data["ElectronAbundance"])
                    T = AU.GetTemp(u, nelec, gamma = 5.0/3.0)
                    rho = AU.PhysicalDensity(np.array(grp_data['Density']),redshift)
                    vol = AU.PhysicalVolume(grp_data['Volume'],redshift)

                    grp_Rvir = grp_data['grp_Rvir']
                    r /= grp_Rvir
                    rho /= AU.ProtonMass

                    for k in np.arange(n_Rbins):
                        in_Rbin = np.logical_and(r > Rbins_min[k], r < Rbins_max[k])
                        z_profile_gal[j,k] = np.sum(m[in_Rbin]*z[in_Rbin])/np.sum(m[in_Rbin])
                        T_profile_gal[j,k] = np.sum(m[in_Rbin]*T[in_Rbin])/np.sum(m[in_Rbin])
                        rho_profile_gal[j,k] = np.sum(m[in_Rbin]*rho[in_Rbin])/np.sum(m[in_Rbin])

                np.savez(self.npz_base+npz_fname,mmin=mmin,mmax=mmax,z_profile_gal=z_profile_gal,T_profile_gal=T_profile_gal,rho_profile_gal=rho_profile_gal)

            if data_type == 'z':
                [Q1,med,Q3] = AU._calc_percentiles(z_profile_gal)
                plt.ylabel(r'$\log_{10}\left[\frac{Z}{Z_\odot}\right]$')
                med /= AU.SolarMetallicity
            if data_type == 'T':
                [Q1,med,Q3] = AU._calc_percentiles(T_profile_gal)
                plt.ylabel(r'$\log_{10}\left[\frac{{\rm Temperature}}{\rm K}\right]$')
            if data_type == 'rho':
                [Q1,med,Q3] = AU._calc_percentiles(rho_profile_gal)
                plt.ylabel(r'Density [cm$^{-3}$]')

            plt.plot(Rbins_min,np.log10(med),label=self.label_list[i],color=self.color_list[i],linestyle=self.linestyle_list[i])

        if data_type == 'z':
            plt.legend(prop={'size':6})
        if data_type == 'T':
            plt.legend(prop={'size':5.5},loc=4,ncol=2,columnspacing=0.3,borderpad=0.2)

        plt.xlabel(r'3D Radius [R$_{200}$]')
        plt.subplots_adjust(left=0.25,bottom=0.18)

        if savename != None:
            plt.savefig(self.fig_base+savename+".pdf")
        else:
            plt.savefig(self.fig_base+"{}_radprof.pdf".format(data_type))




    ##############################
    # 2D data analysis functions #
    ##############################

    def plot_grids(self,species,vmin=10.,vmax=25.,HI_custom_map=False):

        if HI_custom_map:
            # Define custom colormap
            vmin = 10.
            vmax = 25.
            cut_LLS=17
            cut_DLA=20.3
            v_LLS = (cut_LLS-vmin)/(vmax-vmin)
            v_DLA = (cut_DLA-vmin)/(vmax-vmin)

            cdict = {
            'red'  :  ((0., 0., 0.), (v_LLS, 0.3, 1.0),(v_DLA, 0.3, 0.3), (1., 1., 1.)),
            'green':  ((0., 0., 0.), (v_LLS, 0.0, 1.0), (v_DLA, 0.3, 0.0), (1., 0., 0.)),
            'blue' :  ((0., 0., 0.), (v_LLS, 1.0, 0.0), (v_DLA, 0.0, 0.0), (1., 0., 0.))
            }
            cdict = {
            'red'  :  ((0., 0., 0.), (v_LLS, 0.3, 0.3),(v_DLA, 1.0, 0.3), (1., 1., 1.)),
            'green':  ((0., 0., 0.), (v_LLS, 0.0, 0.3), (v_DLA, 1.0, 0.0), (1., 0., 0.)),
            'blue' :  ((0., 0., 0.), (v_LLS, 1.0, 0.0), (v_DLA, 0.0, 0.0), (1., 0., 0.))
            }
            #generate the colormap with 1024 interpolated values
            my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)


        font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 16,
        }
        plt.rc('font',**font)

        # for fn in glob.glob(self.grid_base+"s{}/*.hdf5".format(self.snapnum)):
        for i in xrange(len(self.run_list)):
            run = self.run_list[i]
            snapnum = self.snapnum_list[i]
            # for fn in glob.glob(self.grid_base+"{}/s{}/*.hdf5".format(run,snapnum)):
            for fn in [self.grid_base+"{}/s{}/00032.hdf5".format(run,snapnum)]:
                print "fn ",fn
                f = h5py.File(fn,'r')
                grp_id = f['Header'].attrs['grp_id']
                grid_rad = f['Header'].attrs['grid_radius_pkpc']
                ngrid = f['Header'].attrs['ngrid']
                grid = np.array(f['grids'][species])

                img_savepath = "{}/grids/{}/s{}/{}_{}_300pkpc.pdf".format(self.fig_base,run,snapnum,str(int(grp_id)).zfill(5),species)

                maxdist = grid_rad

                print "grid_rad ",grid_rad
                print "ngrid ",ngrid

                # Check if there are any nans:
                print "Grid has nans? ",np.isnan(np.sum(grid))
                grid[np.isnan(grid)] = 0.
                print "np.min(grid) ",np.min(grid)
                print "np.max(grid) ",np.max(grid)

                plt.close('all')
                if HI_custom_map:
                    plt.imshow(grid,origin='lower',extent=(-maxdist,maxdist,-maxdist,maxdist),vmin=vmin,vmax=vmax,cmap=my_cmap) 
                else: 
                    plt.imshow(grid,origin='lower',extent=(-maxdist,maxdist,-maxdist,maxdist),vmin=vmin,vmax=vmax,cmap=plt.cm.cubehelix) 
                bar=plt.colorbar()
                if species == 'H1':
                    bar_label = r"log$_{10}$ N$_\mathrm{HI}$ (cm$^{-2}$)"
                elif species == 'C4':
                    bar_label = r"log$_{10}$ N$_\mathrm{CIV}$ (cm$^{-2}$)"
                else:
                    bar_label = r"log$_{10}$ N$_\mathrm{"+species+"}$ (cm$^{-2}$)"
                bar.set_label(bar_label)
                plt.xlabel(r"y [proper kpc]")
                plt.ylabel(r"z [proper kpc]")
                plt.savefig(img_savepath, bbox_inches='tight')
                f.close()


    def grid_sightlines(self,species,max_R,coldens_min=0,coldens_max=1000,minmass=10**11.9,maxmass=10**12.1,savename=None,coldens_vs_R=False,coverfrac_vs_R=False,coverfrac_within_R=False,rudie_155=False,rudie_172=False,rudie_19=False,rudie_203=False,turner_c4=False,turner_si4=False,show_Fumagalli=False):
        # Calculates covering fraction in similar way to observers.  Gather all halos within specified mass range, and treat 
        # all corresponding sightlines together.  

        plt.close('all')
        plt.figure(figsize=(3.54331,3.14))

        # Set up bins for radius:
        n_Rbins = 100
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
                # print list(grp_data.keys())
                grid = grp_data[species]
                # Fix for nan "holes":
                grid[np.isnan(grid)] = -1.
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
                [Q1,med,Q3] = AU._calc_percentiles_v2(r,Rbins_min,Rbins_max,N,min_percentile=32,max_percentile=68)

                x = Rbins_med
                x[0] = 0
                x[-1] = max_R
                print "x ",x
                plt.plot(x,med,color=self.color_list[i],label=self.label_list[i],zorder=1)
                if self.run_list[i] == 'gam_50_fixv_nothermal':
                    plt.fill_between(Rbins_med,Q1,Q3,color=self.color_list[i],alpha=0.3,zorder=2)

                if True:
                    thresh = 13.09
                    R_cut = np.logical_and(r>=180.,r<=250.)
                    N_sample = N[R_cut]
                    # now subsample
                    ss = np.zeros(1000)
                    print "subsampling..."
                    for ss_i in np.arange(1000):
                        foo = np.random.choice(N_sample,size=11)
                        ss[ss_i] = np.log10(np.median(10.**foo))
                    print "ss ",ss
                    prob = np.sum(ss >= 13.09)/np.float(np.size(ss))
                    print "Probability of getting observed CIV value at 180-250 kpc: ", prob


            elif coverfrac_vs_R:
                fc = np.zeros(n_Rbins)
                for k in np.arange(n_Rbins):
                    in_Rbin = np.logical_and(r<Rbins_max[k],r>Rbins_min[k])
                    covered = np.logical_and(N[in_Rbin]>=coldens_min,N[in_Rbin]<coldens_max)
                    fc[k] = np.float(np.sum(covered))/np.float(np.sum(in_Rbin))

                plt.plot(Rbins_med,fc,color=self.color_list[i],label=self.label_list[i],zorder=1)
                
            elif coverfrac_within_R:
                fc = np.zeros(n_Rbins)
                for k in np.arange(n_Rbins):
                    in_Rbin = r<Rbins_max[k]
                    covered = np.logical_and(N[in_Rbin]>=coldens_min,N[in_Rbin]<coldens_max)
                    fc[k] = np.float(np.sum(covered))/np.float(np.sum(in_Rbin))

                plt.plot(Rbins_med,fc,color=self.color_list[i],label=self.label_list[i],zorder=1)


        # Rudie et al. 2012, cumulative covering fraction for different values of N_HI for M_halo ~ 10^12
        if rudie_155:
            xdat = np.array([90.,180.])
            ydat = np.array([0.9,0.6])
            yerr = np.array([0.09,0.1])
            xerr = np.array([8.,16.])
            plt.ylim([0.,1.6])
            plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,fmt='.',color='purple',zorder=9)
            if show_Fumagalli:
                plt.scatter(np.array([90.,180.]),np.array([0.38,0.22]),marker='D',color='gray',zorder=10)
        if rudie_172:
            xdat = np.array([90.,180.])
            ydat = np.array([0.3,0.28])
            yerr = np.array([0.14,0.09])
            xerr = np.array([8.,16.])
            plt.ylim([0.,1.0])
            plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,fmt='.',color='purple',zorder=9)
            if show_Fumagalli:
                plt.scatter(np.array([90.,180.]),np.array([0.16,0.07]),marker='D',color='gray',zorder=10)
        elif rudie_19:
            xdat = np.array([90.,180.])
            ydat = np.array([0.1,0.08])
            yerr = np.array([0.09,0.05])
            xerr = np.array([8.,16.])
            plt.ylim([0.,0.4])
            plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,fmt='.',color='purple',zorder=9)
            if show_Fumagalli:
                plt.scatter(np.array([90.,180.]),np.array([0.06,0.03]),marker='D',color='gray',zorder=10)
        elif rudie_203:
            xdat = np.array([90.,180.])
            ydat = np.array([0.,0.04])
            yerr = np.array([0.1,0.04])
            xerr = np.array([8.,16.])
            plt.ylim([0.,0.4])
            plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,fmt='.',color='purple',zorder=9)
            if show_Fumagalli:
                plt.scatter(np.array([90.,180.]),np.array([0.03,0.01]),marker='D',color='gray',zorder=10)


        #Turner+14, column densities vs R for SiIV and CIV:
        if turner_si4:
            pass
        elif turner_c4:
            # impact parameter: 85+-45,155+-25,215+-35,305+-55 
            # column densities:
            x = np.array([85,155,215,305])
            xerr = np.array([45,25,35,55])
            N = np.array([13.62,13.16,13.09,12.66])
            N_ul = np.array([0.12,0.28,0.18,0.36])
            N_ll = np.array([0.18,1.38,0.29,0.5])
            lims = np.array([0,0,0,1])
            # N = np.array([13.55,12.66,12.66,12.80])
            # N_ul = np.array([0.11,0.46,0.30,0.29])
            # N_ll = np.array([0.17,0.5,0.97,1.37])
            # lims = np.array([0,1,0,0])
            plt.errorbar(x,N,xerr=xerr,yerr=[N_ll,N_ul],lolims=lims,linestyle='None',color='purple',zorder=9)
            
        if coldens_vs_R:
            # plt.xlim([0.,max_R])
            plt.xlim([0.,320.])
            plt.ylim([5.,16.])
            print "new ylim!"
            plt.xlabel('Radius [pkpc]')
            if species == 'C4':
                plt.ylabel(r"log$_{10}$ N$_\mathrm{CIV}$ (cm$^{-2}$)")
            elif species == 'O6':
                plt.ylabel(r"log$_{10}$ N$_\mathrm{OVI}$ (cm$^{-2}$)")
            else:
                plt.ylabel(r"log$_{10}$ N$_\mathrm{"+species+"}$ (cm$^{-2}$)")
            # plt.legend(prop={'size':5},ncol=2)
            # plt.ylabel(r"log$_{10}$ N$_\mathrm{OVI}$ (cm$^{-2}$)")
            # plt.xticks(np.array([0.,50.,100.,150.]))
            plt.xticks(np.array([0.,100.,200.,300.]))
            lg=plt.legend(prop={'size':6},loc=3,handletextpad=0.1)
            lg.draw_frame(False)
            plt.title(r"$z=2.5$")
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



    def columndens_vs_R(self,species,kpc_thresh,minmass=0.,maxmass=1e18,plot_mode=1,run_mode='save',savename=None):
        # Plots column density vs radius for a given range of halo masses.  

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
            # foo = "{}kpc".format(kpc_thresh)
            # npz_fname = "fc_{}_s{}_{}_N{}_within{}.npz".format(self.run_list[i],self.snapnum_list[i],species,Nmin,foo)
            # print self.npz_base+npz_fname
            if False: #os.path.isfile(self.npz_base+npz_fname):
                print "Loaded from npz!"
                f = np.load(self.npz_base+npz_fname)
                m = f['m']
                fc_vs_R = f['fc_vs_R']
                f.close()
            else:
                grp_ids = self.find_desired_grpids(i,minmass=minmass,maxmass=maxmass)
                print "grp_ids: ",grp_ids
                print "np.size(grp_ids) ",np.size(grp_ids)
                m = np.zeros(np.size(grp_ids))
                fc = np.zeros(np.size(grp_ids))
                N_vs_R = np.zeros([np.size(grp_ids),n_Rbins])
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

                    if kpc_thresh > grid_rad:
                        raise Exception('Desired radial threshold ({}) exceeded grid radius ({}) for covering fraction calculation.'.format(dist_thresh,grid_rad))

                    for R_ind in np.arange(n_Rbins):
                        in_Rbin = r_kpc<=Rbins_max[R_ind]
                        if np.sum(in_Rbin) == 0:
                            N_vs_R[j,R_ind] = 0.
                        else:
                            #in_N_range = np.logical_and(grid[in_Rbin]>Nmin,grid[in_Rbin]<Nmax)
                            #fc_vs_R[j,R_ind] = np.float(np.sum(in_N_range))/np.float(np.sum(in_Rbin))
                            N_vs_R[j,R_ind] = np.mean(grid[in_Rbin]) # effectively averaging by area instead of by mass

                # Save covering fraction data to npz file
                # np.savez(self.npz_base+npz_fname,run=self.run_list[i],snapnum=self.snapnum_list[i],species=species,Nmin=Nmax,kpc_thresh=kpc_thresh,grp_ids=grp_ids,m=m,fc_vs_R=fc_vs_R)  

            [Q1,med,Q3] = AU._calc_percentiles(N_vs_R,min_percentile=10,max_percentile=90)

            if plot_mode == 1:
                plt.plot(Rbins_min,med,color=self.color_list[i],label=self.label_list[i],linestyle=self.linestyle_list[i])
                if plot_mode == 2:
                    plt.fill_between(Rbins_min,Q1,Q3,color=self.color_list[i],alpha=0.3)
                plt.xlabel('Radius [pkpc]')
                plt.ylabel(r'{} $N(< r)$'.format(species))
                plt.legend()

        if savename != None:
            plt.savefig(self.fig_base+savename+".pdf")
        else:
            plt.savefig(self.fig_base+"N_vs_R.pdf")
                

    def coverfrac_bootstrap(self,species,max_R,coldens_min=0,coldens_max=1000,minmass=10**11.9,maxmass=10**12.1,savename=None,rudie_155=False,rudie_172=False,rudie_19=False,rudie_203=False,show_Fumagalli=False):
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
                # print list(grp_data.keys())
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


        # Bootstrap for < 90 kpc, and 10 sightlines:
        print "Starting bootstrap now..."
        N_sample = 100000
        r_cut = r <= 180.
        N_dist = N[r_cut]
        fc_dist = np.zeros(N_sample)
        for i in np.arange(N_sample):
            s = np.random.choice(N_dist,size=25)
            fc_dist[i] = np.sum(np.logical_and(s > coldens_min, s < coldens_max))/25.
        # Now we find the median and quartiles of distribution:
        print "fc_dist ",fc_dist
        print "Bootstrap results: "
        print "Median: ",np.median(fc_dist)
        print "1-sig below: ", np.percentile(fc_dist,32)
        print "1-sig above: ", np.percentile(fc_dist,68)
        print "2-sig below: ", np.percentile(fc_dist,5)
        print "2-sig above: ", np.percentile(fc_dist,95)
        print "std-dev: ", np.std(fc_dist)
        print "mean: ",np.mean(fc_dist)

        # Bootstrap for < 180 kpc, and 25 sightlines.  10 within 90 kpc, 15 within 90-180 kpc:
        print "Starting bootstrap now..."
        N_sample = 100000
        r_cut1 = r <= 90.
        r_cut2 = np.logical_and(r > 90.,r <= 180.)
        N_dist1 = N[r_cut1]
        N_dist2 = N[r_cut2]
        fc_dist = np.zeros(N_sample)
        for i in np.arange(N_sample):
            s = np.append(np.random.choice(N_dist1,size=10),np.random.choice(N_dist2,size=15))
            fc_dist[i] = np.sum(np.logical_and(s > coldens_min, s < coldens_max))/25.
        # Now we find the median and quartiles of distribution:
        print "fc_dist ",fc_dist
        print "Bootstrap results: "
        print "Median: ",np.median(fc_dist)
        print "1-sig below: ", np.percentile(fc_dist,32)
        print "1-sig above: ", np.percentile(fc_dist,68)
        print "2-sig below: ", np.percentile(fc_dist,5)
        print "2-sig above: ", np.percentile(fc_dist,95)
        print "mean ",np.mean(fc_dist)
        print "std-dev: ", np.std(fc_dist)



        if False:
            plt.plot(Rbins_med,fc,color=self.color_list[i],label=self.label_list[i],zorder=1)


            # Rudie et al. 2012, cumulative covering fraction for different values of N_HI for M_halo ~ 10^12
            if rudie_155:
                xdat = np.array([90.,180.])
                ydat = np.array([0.9,0.6])
                yerr = np.array([0.09,0.1])
                xerr = np.array([8.,16.])
                plt.ylim([0.,1.6])
                plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,fmt='.',color='purple',zorder=9)
                if show_Fumagalli:
                    plt.scatter(np.array([90.,180.]),np.array([0.38,0.22]),marker='D',color='gray',zorder=10)
            if rudie_172:
                xdat = np.array([90.,180.])
                ydat = np.array([0.3,0.28])
                yerr = np.array([0.14,0.09])
                xerr = np.array([8.,16.])
                plt.ylim([0.,1.0])
                plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,fmt='.',color='purple',zorder=9)
                if show_Fumagalli:
                    plt.scatter(np.array([90.,180.]),np.array([0.16,0.07]),marker='D',color='gray',zorder=10)
            elif rudie_19:
                xdat = np.array([90.,180.])
                ydat = np.array([0.1,0.08])
                yerr = np.array([0.09,0.05])
                xerr = np.array([8.,16.])
                plt.ylim([0.,0.4])
                plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,fmt='.',color='purple',zorder=9)
                if show_Fumagalli:
                    plt.scatter(np.array([90.,180.]),np.array([0.06,0.03]),marker='D',color='gray',zorder=10)
            elif rudie_203:
                xdat = np.array([90.,180.])
                ydat = np.array([0.,0.04])
                yerr = np.array([0.1,0.04])
                xerr = np.array([8.,16.])
                plt.ylim([0.,0.4])
                plt.errorbar(xdat,ydat,yerr=yerr,xerr=xerr,fmt='.',color='purple',zorder=9)
                if show_Fumagalli:
                    plt.scatter(np.array([90.,180.]),np.array([0.03,0.01]),marker='D',color='gray',zorder=10)


                
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

    def find_desired_grpids(self,i,minmass=0,maxmass=1e20,return_grpmass=False):
        print self.snapdir_list[i],self.snapnum_list[i]
        cat = readsubfHDF5.subfind_catalog(self.snapdir_list[i],self.snapnum_list[i],subcat=False)
        grp_mass = AU.PhysicalMass(np.array(cat.Group_M_Crit200))
        mass_select = np.logical_and(grp_mass > minmass,grp_mass < maxmass)

        grp_ids = np.arange(np.float64(cat.ngroups))
        grp_ids = grp_ids[mass_select]
        # grp_mass = grp_mass[mass_select]
        # for j in np.arange(np.size(grp_ids)):
        #     print "grp id {} - grp mass {}".format(grp_ids[j],grp_mass[j])
        if not return_grpmass:
            return np.sort(grp_ids)
        else:
            foo = np.argsort(grp_ids)
            return [grp_ids[foo],grp_mass[foo]]


    def get_gal_props(self,run_index,grp_ids):
        snapdir = self.snapdir_list[run_index]
        snapnum = self.snapnum_list[run_index]
        cat = readsubfHDF5.subfind_catalog(snapdir,snapnum,long_ids=True)

        grp_ids = np.int32(grp_ids)
        sub_ids = cat.GroupFirstSub[grp_ids]
        gal_mass = AU.PhysicalMass(cat.SubhaloMassType[sub_ids,4])
        # gal_SFR = AU.PhysicalMdot(cat.SubhaloSFR[sub_ids])
        # data_dict = {'gal_mass':gal_mass,'gal_SFR':gal_SFR}
        data_dict = {'gal_mass':gal_mass}
        return data_dict

    def load_CGMsnap_data(self,run_index,grp_id,data_entry='all'):
        fname = "{}{}/s{}/{}.hdf5".format(self.CGMsnap_base,self.run_list[run_index],self.snapnum_list[run_index],str(int(grp_id)).zfill(5))
        #fname = self.CGMsnap_base + self.run_list[run_index] + str(grp_id).zfill(5) + ".hdf5"
        # check that file exists - NOT IMPLEMENTED
        print "fname: ",fname
        f = h5py.File(fname,'r')

        # if data_entry is 'all', return a dictionary with all data fields, otherwise return specific data field
        if data_entry == 'all':
            data_dict = {}
            for key in f['Header'].attrs:
                data_dict[key] = f['Header'].attrs[key]
            for key in f['PartType0']:
                data_dict[key] = f['PartType0'][key]
            return data_dict
        else:
            return f['PartType0'][data_entry]

    def load_grid_data(self,i,grp_id):
        # assume that i is singular
        grp_data = {}
        # print self.grid_base+"{}/s{}/{}.hdf5".format(self.run_list[i],self.snapnum_list[i],str(int(grp_id)).zfill(5))
        try: f = h5py.File(self.grid_base+"{}/s{}/{}.hdf5".format(self.run_list[i],self.snapnum_list[i],str(int(grp_id)).zfill(5)),'r')
        except: 
            print "File could not be opened!: {}".format(self.grid_base+"{}/s{}/{}.hdf5".format(self.run_list[i],self.snapnum_list[i],str(int(grp_id)).zfill(5)))
        for key in list(f['Header'].attrs): grp_data[key] = f['Header'].attrs[key]
        for key in f['grids'].keys(): grp_data[key] = np.array(f['grids'][key])
        print "returning grp data from ",self.grid_base+"{}/s{}/{}.hdf5".format(self.run_list[i],self.snapnum_list[i],str(int(grp_id)).zfill(5))
        return grp_data 

    def _grid_to_kpc(self,r_grid,ngrid,grid_rad):
        return r_grid * (2.*grid_rad)/(ngrid)

    def _calc_radii(self,grp_data):
        box = grp_data['box']
        grp_pos = grp_data['grp_pos']
        pos = grp_data['Coordinates']
        r = AU._pbc_dist(pos,grp_pos,boxsize=box)
        return r
        
    # self.get_gal_props


    def dat_prep(self,redshifts,c0_256=0,c0_512=0,c0_128=0,c3_512=0,c4_512=0,c4_check=0,c2_256=0,c0_fw_256=0,c0_sw_256=0,c5_256=0,g0_BH=0,g10_BH=0,g20_BH=0,g25_BH=0,g30_BH=0,g40_BH=0,g25_noBH=0,g50_BH=0,g75_BH=0,g75_noBH=0,g95_BH=0,g95_noBH=0,g10_noBH=0,g20_noBH=0,g30_noBH=0,g40_noBH=0,g50_noBH=0,g10_nothermal=0,g20_nothermal=0,g30_nothermal=0,g40_nothermal=0,g50_nothermal=0,g25_fixv=0,g50_fixv=0,g25_fixv_nothermal=0,g50_fixv_nothermal=0,g25_fixv_fixeta=0,g50_fixv_fixeta=0,c0_nometalwinds=0,c0_fullmetalwinds=0,c0_dumpfactor95=0,no_metal_cooling=0,low_recouple_dens=0,high_recouple_dens=0,multi_speed=0):
        mv_snapbase = "/n/hernquistfs1/Illustris/SmallBox/GFM/Production/Cosmo/" 
        sb_snapbase = "/n/hernquistfs1/spb/Cosmo/"
        js_snapbase = "/n/home04/jsuresh/runs/Cosmo5_V6/output"
        # gam_snapbase = "/n/hernquistfs1/jsuresh/Runs/"
        gam_snapbase = "/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/Runs/"
        ap_snapbase = "/n/regal/hernquist_lab/apillepich/Simulations_2014_IllustrisTNG_Runs/"

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
                    self.label_list.append("Fiducial")
                    self.snapnum_list.append(60)
                elif redshift == '2.5':
                    self.label_list.append("Fiducial")
                    self.snapnum_list.append(63)
                elif redshift == '2':
                    self.label_list.append(r"$256^3$")
                    # self.label_list.append("Fiducial")
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
                #self.color_list.append(cmap(0.))
                self.color_list.append('red')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(mv_snapbase+"Cosmo0_V6/L25n512/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0+BHs (z=4)")
                    self.snapnum_list.append(54)
                elif redshift == '3':
                    self.label_list.append("Fiducial")
                    self.snapnum_list.append(60)
                elif redshift == '2.5':
                    self.label_list.append("Gamma=0+BHs (z=2.5)")
                    self.snapnum_list.append(63)
                elif redshift == '2':
                    self.label_list.append(r"$512^3$")
                    # self.label_list.append("Fiducial")
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

        if c0_128:
            for redshift in redshifts:
                self.run_list.append("c0_128")
                #self.color_list.append(cmap(0.))
                self.color_list.append('green')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(mv_snapbase+"Cosmo0_V6/L25n128/output/")
                if redshift == '4':
                    self.label_list.append("Gamma=0+BHs (z=4)")
                    self.snapnum_list.append(54)
                elif redshift == '3':
                    self.label_list.append("Fiducial")
                    self.snapnum_list.append(60)
                elif redshift == '2.5':
                    self.label_list.append("Gamma=0+BHs (z=2.5)")
                    self.snapnum_list.append(63)
                elif redshift == '2':
                    self.label_list.append(r"$128^3$")
                    # self.label_list.append("Fiducial")
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
                    self.label_list.append("No AGN")
                    self.snapnum_list.append(60)
                elif redshift == '2':
                    # self.label_list.append("Cold Winds only (z=2)")
                    self.label_list.append("No AGN")
                    self.snapnum_list.append(68)
                elif redshift == '2.5':
                    # self.label_list.append("Cold Winds only (z=2)")
                    self.label_list.append("No AGN")
                    self.snapnum_list.append(63)
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
                    self.label_list.append("Higher Mass-loading")
                    self.snapnum_list.append(60)
                elif redshift == '2':
                    self.label_list.append("Higher Mass-loading")
                    self.snapnum_list.append(68)
                elif redshift == '2.5':
                    self.label_list.append("Higher Mass-loading")
                    self.snapnum_list.append(63)
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
                    self.label_list.append("Fiducial")
                    self.snapnum_list.append(60)
                elif redshift == '2.5':
                    self.label_list.append("Fiducial")
                    self.snapnum_list.append(63)
                elif redshift == '2':
                    self.label_list.append(r"$256^3$")
                    # self.label_list.append("Fiducial")
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
                    self.label_list.append("Fixed-E Hot Winds")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Fixed-E Hot Winds")
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
                    self.label_list.append("Fixed-v Hot Winds")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("Fixed-v Hot Winds")
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
                if redshift == '3':
                    self.label_list.append("Faster Winds")
                    self.snapnum_list.append(3)
                if redshift == '2.5':
                    self.label_list.append("Faster Winds")
                    self.snapnum_list.append(4)
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

        if c0_dumpfactor95:
            for redshift in redshifts:
                self.run_list.append("c0_dumpfactor95")
                self.color_list.append("green")
                self.linestyle_list.append("solid")
                self.snapdir_list.append(gam_snapbase+"c0_dumpfactor95/output/")
                if redshift == '4':
                    self.label_list.append("95\% Dump Factor")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("95\% Dump Factor")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("95\% Dump Factor")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    self.label_list.append("95\% Dump Factor")
                    self.snapnum_list.append(5)

        if no_metal_cooling:
            for redshift in redshifts:
                self.run_list.append("no_metal_cooling")
                self.color_list.append("indigo")
                self.linestyle_list.append("solid")
                self.snapdir_list.append(gam_snapbase+"no_metal_cooling/output/")
                if redshift == '4':
                    self.label_list.append("No Metal Cooling")
                    self.snapnum_list.append(1)
                elif redshift == '3':
                    self.label_list.append("No Metal Cooling")
                    self.snapnum_list.append(3)
                elif redshift == '2.5':
                    self.label_list.append("No Metal Cooling")
                    self.snapnum_list.append(4)
                elif redshift == '2':
                    self.label_list.append("No Metal Cooling")
                    self.snapnum_list.append(5)

        if low_recouple_dens:
            for redshift in redshifts:
                self.run_list.append("low_recouple_dens")
                self.color_list.append("magenta")
                self.linestyle_list.append("solid")
                self.snapdir_list.append("/n/regal/hernquist_lab/apillepich/Simulations_2014_IllustrisTNG_Runs/APTests_L25n256FP_0000_2501/Output/")
                if redshift == '2':
                    self.label_list.append("Low Recoupling Density")
                    self.snapnum_list.append(17)

        if high_recouple_dens:
            for redshift in redshifts:
                self.run_list.append("high_recouple_dens")
                self.color_list.append("cyan")
                self.linestyle_list.append("solid")
                self.snapdir_list.append("/n/regal/hernquist_lab/apillepich/Simulations_2014_IllustrisTNG_Runs/APTests_L25n256FP_0000_2502/Output/")
                if redshift == '2':
                    self.label_list.append("High Recoupling Density")
                    self.snapnum_list.append(17)

        if multi_speed:
            for redshift in redshifts:
                self.run_list.append("multi_speed")
                self.color_list.append("red")
                self.linestyle_list.append("solid")
                self.snapdir_list.append("/n/home01/ptorrey/Share/apillepich/L25n256/multi_speed_winds/output/")
                if redshift == '2.5':
                    self.label_list.append("Multi-speed")
                    self.snapnum_list.append(15)
                if redshift == '2':
                    self.label_list.append("Multi-speed")
                    self.snapnum_list.append(17)

        self.ndat = len(self.run_list)

if __name__ == '__main__':
    full_analysis()

