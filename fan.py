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
        # self.obs_cover_frac_vs_R("H1",17.2,200.,Nmax=20.3,minmass=10**11.8,maxmass=10**12.2,savename="LLS_postisaac",rudie_LLS=True)
        # self.obs_cover_frac_vs_R("H1",20.3,200.,minmass=10**11.8,maxmass=10**12.2,savename="DLA_postisaac",rudie_DLA=True)
        # add a column-density vs R function eventually

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

    #############
    # Utilities #
    #############
    # self.find_desired_grpids
    # self.load_grid_data
    # self._grid_to_kpc
    # self.get_gal_props

    def dat_prep(self,redshifts,c0_128=0,c0_256=0,c0_512=0,c0_fw_256=0,c0_sw_256=0,c2_256=0,c3_512=0,c4_512=0,c4_check=0,g0_BH=0,g10_BH=0,g20_BH=0,g30_BH=0,g40_BH=0,g50_BH=0,g25_BH=0,g75_BH=0,g95_BH=0,g0_noBH=0,g10_noBH=0,g20_noBH=0,g30_noBH=0,g40_noBH=0,g50_noBH=0,g25_noBH=0,g75_noBH=0,g95_noBH=0,g10_nothermal=0,g20_nothermal=0,g30_nothermal=0,g40_nothermal=0,g50_nothermal=0,g25_fixv=0,g50_fixv=0,g25_fixv_fixeta=0,g50_fixv_fixeta=0,g25_fixv_nothermal=0,g50_fixv_nothermal=0,c0_nometalwinds=0,c0_fullmetalwinds=0):
        # self.run_list = []
        # self.color_list = [] 
        # self.label_list = [] 
        # self.linestyle_list = []
        # self.snapnum_list = [] 
        # self.snapdir_list = []

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
                self.color_list.append('blue')
                self.label_list.append('Strong Winds')
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
                self.color_list.append('blue')
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
                self.color_list.append('blue')
                self.label_list.append('Fiducial')
                self.linestyle_list.append('solid')
                self.snapdir_list.append(self.mv_base+"Cosmo0_V6/L25_n256")
                if redshift == '4': self.snapnum_list.append(54)
                if redshift == '3': self.snapnum_list.append(60)
                if redshift == '2': self.snapnum_list.append(68)

        # self.dat_prep(redshifts,g0_BH=1)
        # self.dat_prep(redshifts,c2_256=1)
        # self.dat_prep(redshifts,c0_sw_256=1)
        # self.dat_prep(redshifts,g50_fixv_nothermal=1)
        # self.dat_prep(redshifts,g50_BH=1)
        # self.dat_prep(redshifts,g50_fixv=1)
        # self.dat_prep(redshifts,c0_nometalwinds=1)
        # self.dat_prep(redshifts,c0_fullmetalwinds=1)













      