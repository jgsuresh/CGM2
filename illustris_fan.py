import numpy as np
import readsubfHDF5
import h5py

import matplotlib
matplotlib.use('PDF')
from matplotlib import rc
#rc('text', usetex=True)
rc('font', family='serif')
import matplotlib.pyplot as plt
import matplotlib.ticker
import brewer2mpl

from units import AREPO_units
AU = AREPO_units()

import convert_cloudy as cc

import glob
import os

class illustris_fan:
    def __init__(self):
        self.snapbase = "/n/ghernquist/Illustris/Runs/Illustris-1/output/"
        # self.CGMsnap_base = '/n/home04/jsuresh/scratch1/AREPOfest/data/CGM_snaps/'
        self.CGMsnap_base = '/n/home04/jsuresh/scratch1/CGM_snaps/'
        self.grid_base = '/n/home04/jsuresh/scratch1/AREPOfest/data/grids/'
        self.fig_base = '/n/home04/jsuresh/scratch1/AREPOfest/data/figs/'
        self.npz_base = '/n/home04/jsuresh/scratch1/AREPOfest/data/npz/'
        # self.CGMsnap_base = '/n/home04/jsuresh/scratch1/QCGM2/data/CGM_snaps/'
        # self.grid_base = '/n/home04/jsuresh/scratch1/QCGM2/data/grids/'
        # self.fig_base = '/n/home04/jsuresh/scratch1/QCGM2/data/figs/'
        # self.npz_base = '/n/home04/jsuresh/scratch1/QCGM2/data/npz/'

        self.snapnum = 120
        self.redshift = 0.2

        #read in catalog + galprop file
        # self.cat = readsubfHDF5.subfind_catalog('/n/ghernquist/Illustris/Runs/Illustris-1/',self.snapnum,keysel=['GroupFirstSub'],subcat=False)
        self.galf = h5py.File(self.snapbase+"/postprocessing/galprop/galprop_{}.hdf5".format(str(self.snapnum).zfill(3)),'r')

        self.cat = readsubfHDF5.subfind_catalog(self.snapbase, self.snapnum, long_ids=True, double_output=True, keysel=["GroupFirstSub","SubhaloGrNr"])
        self.load_gal_props()
        # self.gal_mass_vs_sSFR()
        # self.halo_baryon_abundance()

        self.tab = cc.CloudyTable(0.2,directory='/n/home04/jsuresh/scratch1/Cloudy_test/UVB_sf_xrays_ext/')
        self.OVI_weighted_gas('T',savename='T_test_v2')
        # self.OVI_weighted_gas('rho',savename='rho_test')

        # self.plot_grids("H1",vmax=25.)
        # self.plot_grids("Si3",vmax=18.)
        # self.plot_grids("N5",vmax=16.)
        # self.plot_grids("O6",vmax=16.)
        # self.galprop_vs_CGMprop('sm','CGM_ISM_metal_ratio')
        # self.gal_mass_vs_Rvir()
        # self.coldens_plot('N5',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,savename='N5_test')
        # self.coldens_plot('N5',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='UVB_sf_xrays_ext',savename='N5_lowssfr_medm_xrays')
        # self.coldens_plot('N5',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='ion_out_fancy_atten',savename='N5_lowssfr_medm_UVB')
        # self.coldens_plot('N5',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='UVB_sf_xrays_ext',savename='N5_highssfr_medm_xrays')
        # self.coldens_plot('N5',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='ion_out_fancy_atten',savename='N5_highssfr_medm_UVB')
        # self.coldens_plot('N5',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='UVB_sf_xrays_ext',savename='N5_lowssfr_highm_xrays')
        # self.coldens_plot('N5',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='ion_out_fancy_atten',savename='N5_lowssfr_highm_UVB')
        # self.coldens_plot('N5',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='UVB_sf_xrays_ext',savename='N5_highssfr_highm_xrays')
        # self.coldens_plot('N5',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='ion_out_fancy_atten',savename='N5_highssfr_highm_UVB')

        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='UVB_sf_xrays_ext',tempfac=1,savename='O6_highssfr_medm_xrays_t1_test')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='UVB_sf_xrays_ext',tempfac=1,savename='O6_lowssfr_medm_xrays_t1_test')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='UVB_sf_xrays_ext',tempfac=1,savename='N5_highssfr_medm_xrays')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='ion_out_fancy_atten',tempfac=1,savename='N5_highssfr_medm_UVB')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='UVB_sf_xrays_ext',tempfac=1,savename='N5_lowssfr_highm_xrays')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='ion_out_fancy_atten',tempfac=1,savename='N5_lowssfr_highm_UVB')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='UVB_sf_xrays_ext',tempfac=1.,savename='O6_highssfr_highm_xrays_t1')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='ion_out_fancy_atten',tempfac=1.,savename='O6_highssfr_highm_UVB_t1')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='UVB_sf_xrays_ext',tempfac=3.,savename='O6_highssfr_highm_xrays_t3')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='ion_out_fancy_atten',tempfac=3.,savename='O6_highssfr_highm_UVB_t3')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='UVB_sf_xrays_ext',tempfac=10.,savename='O6_highssfr_highm_xrays_t10')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='ion_out_fancy_atten',tempfac=10.,savename='O6_highssfr_highm_UVB_t10')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='UVB_sf_xrays_ext',tempfac=1.,savename='O6_lowssfr_highm_xrays_t1')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='ion_out_fancy_atten',tempfac=1.,savename='O6_lowssfr_highm_UVB_t1')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='UVB_sf_xrays_ext',tempfac=3.,savename='O6_lowssfr_highm_xrays_t3')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='ion_out_fancy_atten',tempfac=3.,savename='O6_lowssfr_highm_UVB_t3')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='UVB_sf_xrays_ext',tempfac=10.,savename='O6_lowssfr_highm_xrays_t10')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**11.,Mmax=10.**11.5,cloudy='ion_out_fancy_atten',tempfac=10.,savename='O6_lowssfr_highm_UVB_t10')

        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='UVB_sf_xrays_ext',tempfac=1.,savename='O6_lowssfr_medm_xrays_t1')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='ion_out_fancy_atten',tempfac=1.,savename='O6_lowssfr_medm_UVB_t1')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='UVB_sf_xrays_ext',tempfac=3.,savename='O6_lowssfr_medm_xrays_t3')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='ion_out_fancy_atten',tempfac=3.,savename='O6_lowssfr_medm_UVB_t3')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='UVB_sf_xrays_ext',tempfac=10.,savename='O6_lowssfr_medm_xrays_t10')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='ion_out_fancy_atten',tempfac=10.,savename='O6_lowssfr_medm_UVB_t10')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='UVB_sf_xrays_ext',tempfac=3.,savename='O6_highssfr_medm_xrays_t3')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='ion_out_fancy_atten',tempfac=3.,savename='O6_highssfr_medm_UVB_t3')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='UVB_sf_xrays_ext',tempfac=10.,savename='O6_highssfr_medm_xrays_t10')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='ion_out_fancy_atten',tempfac=10.,savename='O6_highssfr_medm_UVB_t10')

        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='UVB_sf_xrays_ext',savename='O6_lowssfr_medm_xrays_mult')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,low_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='ion_out_fancy_atten',savename='O6_lowssfr_medm_UVB_mult')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='UVB_sf_xrays_ext',savename='O6_highssfr_medm_xrays_mult')
        # self.coldens_plot('O6',vmin=11.,vmax=18,kpc_mode=True,high_ssfr_pop=True,Mmin=10.**10.5,Mmax=10.**11.,cloudy='ion_out_fancy_atten',savename='O6_highssfr_medm_UVB_mult')

        # self.mass_metal_budget(BH_split=True)

        self.galf.close()


    def gal_mass_vs_sSFR(self,include_COS=True,savename=None):
        plt.figure()

        ax = plt.scatter(np.log10(self.gal_sm),np.log10(self.gal_ssfr),marker='.',s=10,alpha=0.3,zorder=1)

        if include_COS:
            werk_dat = np.loadtxt("werk.dat")
            COS_log_sm = werk_dat[:,0]
            COS_sfr = werk_dat[:,1]
            upper_lim = np.array(werk_dat[:,2],dtype='bool')
            not_upper_lim = np.logical_not(upper_lim)

            COS_sm = 10.**COS_log_sm
            COS_ssfr = COS_sfr/COS_sm
            plt.scatter(np.log10(COS_sm[not_upper_lim]),np.log10(COS_ssfr[not_upper_lim]),marker='o',color='purple',zorder=2)
            plt.scatter(np.log10(COS_sm[upper_lim]),np.log10(COS_ssfr[upper_lim]),marker='v',color='purple',zorder=3)

        # ax.set_xscale('log')
        # ax.set_yscale('log')
        if savename != None: 
            plt.savefig(self.fig_base+savename+".pdf", bbox_inches='tight')
        else:
            plt.savefig(self.fig_base+"sm_vs_ssfr.pdf", bbox_inches='tight')


    def gal_mass_vs_Rvir(self,savename=None):
        plt.figure()

        self.load_CGMsnap_ids()
        print "about to get cat "
        cat = readsubfHDF5.subfind_catalog('/n/ghernquist/Illustris/Runs/Illustris-1/',self.snapnum,keysel=['Group_M_Crit200','Group_R_Crit200','GroupFirstSub'],subcat=False)
        print "got cat"

        self.sub_ids = cat.GroupFirstSub[np.int32(self.grp_ids)]
        galf = h5py.File(self.snapbase+"/postprocessing/galprop/galprop_{}.hdf5".format(self.snapnum),'r')
        gal_sm = AU.PhysicalMass(np.array(galf['stellar_totmass']))
        x = gal_sm[self.sub_ids]

        Rvir = AU.PhysicalPosition(cat.Group_M_Crit200[np.int32(self.grp_ids)],0.19728)
        # plt.scatter(np.log10(x),Rvir,marker='.',s=10,alpha=0.3,zorder=1)
        plt.scatter(np.log10(x),np.log10(AU.PhysicalMass(cat.Group_M_Crit200[np.int32(self.grp_ids)])),marker='.',s=10,alpha=0.3,zorder=1)
        # plt.vline(150.)

        # ax.set_xscale('log')
        # ax.set_yscale('log')
        if savename != None: 
            plt.savefig(self.fig_base+savename+".pdf", bbox_inches='tight')
        else:
            plt.savefig(self.fig_base+"sm_vs_halomass.pdf", bbox_inches='tight')


    def halo_baryon_abundance(self,savename=None):
        self.load_CGMsnap_ids()

        print "self.grp_ids ",self.grp_ids
        cat = readsubfHDF5.subfind_catalog('/n/ghernquist/Illustris/Runs/Illustris-1/',self.snapnum,keysel=['Group_M_Crit200'],subcat=False)
        grp_mass = cat.Group_M_Crit200[np.int32(self.grp_ids)]

        fb = (0.0456/0.27)
        gas_mass = np.zeros(0)
        for grp_id in self.grp_ids:
            f = self.CGMsnap_base+"s{}/{}.hdf5".format(self.snapnum,str(int(grp_id)).zfill(5))
            data_dict = self.load_CGM_snap(f)
            gas_mass = np.append(gas_mass,np.sum(data_dict['MASS']))

        frac = gas_mass/(fb*grp_mass)

        plt.figure()
        x = np.log10(AU.PhysicalMass(grp_mass))
        y = frac
        plt.scatter(x,y)
        plt.savefig(self.fig_base+"halo_gas_abundance.pdf")

        plt.figure()
        x = np.log10(AU.PhysicalMass(grp_mass))
        y = gas_mass/grp_mass
        plt.scatter(x,y)
        plt.xlabel(r"$M_{200}$")
        plt.ylabel(r"$\frac{M_{\rm gas}}{M_{200}}$")
        plt.savefig(self.fig_base+"halo_gas_abundance_v2.pdf")


    def mass_metal_budget(self,savename=None,BH_split=False):

        # # Get CGMsnap group ids, and main subhalo ids:
        # self.load_CGMsnap_ids()
        # cat = readsubfHDF5.subfind_catalog('/n/ghernquist/Illustris/Runs/Illustris-1/',self.snapnum,keysel=['GroupPos','GroupFirstSub','Group_M_Crit200','SubhaloStarMetallicity'])
        # self.sub_ids = np.int32(cat.GroupFirstSub[np.int32(self.grp_ids)])

        f_npz = self.npz_base + "{}.npz".format("mass_met_budget_Rvir")
        if os.path.isfile(f_npz):
            f_npz = np.load(f_npz)
            m_grp = f_npz['m_grp']
            m_winds = f_npz['m_winds']
            mz_winds = f_npz['mz_winds']
            m_stars = f_npz['m_stars']
            mz_stars = f_npz['mz_stars']
            m_ISM = f_npz['m_ISM']
            mz_ISM = f_npz['mz_ISM']
            m_CGM = f_npz['m_CGM']
            mz_CGM = f_npz['mz_CGM']
            m_cool_CGM = f_npz['m_cool_CGM']
            mz_cool_CGM = f_npz['mz_cool_CGM']
            m_warm_CGM = f_npz['m_warm_CGM']
            mz_warm_CGM = f_npz['mz_warm_CGM']
            m_hot_CGM = f_npz['m_hot_CGM']
            mz_hot_CGM = f_npz['mz_hot_CGM']
        else:
            # Get CGMsnap group ids, and main subhalo ids:
            self.load_CGMsnap_ids()
            cat = readsubfHDF5.subfind_catalog('/n/ghernquist/Illustris/Runs/Illustris-1/',self.snapnum,keysel=['GroupPos','GroupFirstSub','Group_M_Crit200','SubhaloStarMetallicity'])
            self.sub_ids = np.int32(cat.GroupFirstSub[np.int32(self.grp_ids)])

            
            # # check Wind/Stars mass ratio
            # gal_wm = AU.PhysicalMass(np.array(self.galf['wind_totmass']))
            # gal_sm = AU.PhysicalMass(np.array(self.galf['stellar_totmass']))
            # print "gal_wm/gal_sm ",gal_wm/gal_sm 
            # # i = gal_sm > 0
            # f = np.log10(gal_wm[self.sub_ids]/gal_sm[self.sub_ids])
            # f = f[f==f]
            # print "np.median(np.log10(gal_wm/gal_sm)) ",np.median(f)

            # Group:
            m_grp = AU.PhysicalMass(cat.Group_M_Crit200[np.int32(self.grp_ids)])

            # Winds:
            m_winds = AU.PhysicalMass(np.array(self.galf['wind_totmass'])[self.sub_ids])
            mz_winds = 0.4*cat.SubhaloStarMetallicity[self.sub_ids]*m_winds # VERY ROUGH ESTIMATE!

            # Stars:
            m_stars = AU.PhysicalMass(np.array(self.galf['stellar_totmass'])[self.sub_ids])
            mz_stars = cat.SubhaloStarMetallicity[self.sub_ids]*m_stars
            gal_ssfr = np.array(self.galf['gas_totsfr'])[self.sub_ids]/m_stars

            # ISM and CGM
            def load_CGM_snap(fn,load='all'):
                data_dict = {}
                f = h5py.File(fn,'r')
                if load=='all':
                    for key in f['Header'].attrs:
                        data_dict[key] = f['Header'].attrs[key]
                    for key in f['PartType0']:
                        data_dict[key] = np.copy(f['PartType0'][key])
                else:
                    data_dict[key] = f['Header'].attrs[key]
                return data_dict

            i = 0
            m_ISM = np.zeros(0) # from CGM_snap, with radial cut for 150 kpc
            mz_ISM = np.zeros(0) # from CGM_snap, with radial cut for 150 kpc
            m_CGM = np.zeros(0) # from CGM_snap, with radial cut for 150 kpc
            mz_CGM = np.zeros(0) # from CGM_snap, with radial cut for 150 kpc
            m_cool_CGM = np.zeros(0) # from CGM_snap, with radial cut for 150 kpc and T < 10^5 K.
            mz_cool_CGM = np.zeros(0) # from CGM_snap, with radial cut for 150 kpc and T < 10^5 K.
            m_warm_CGM = np.zeros(0) # from CGM_snap, with radial cut for 150 kpc and 10^5 K < T < 10^6 K.
            mz_warm_CGM = np.zeros(0) # from CGM_snap, with radial cut for 150 kpc and 10^5 K < T < 10^6 K.
            m_hot_CGM = np.zeros(0) # from CGM_snap, with radial cut for 150 kpc and T > 10^6 K.
            mz_hot_CGM = np.zeros(0) # from CGM_snap, with radial cut for 150 kpc and T > 10^6 K.

            # for fn in glob.glob(self.CGMsnap_base+"s{}/*.hdf5".format(self.snapnum)):
            # self.grp_ids = self.grp_ids[:25]
            for grp_id in self.grp_ids:
                fn = self.CGMsnap_base+"s{}/{}.hdf5".format(self.snapnum,str(int(grp_id)).zfill(5))
                print "fn ",fn
                print "i ",i
                data_dict = load_CGM_snap(fn)
                rho = AU.PhysicalDensity(np.array(data_dict["RHO "]),0.19728)
                x = np.array(data_dict["POS "])

                grp_pos = cat.GroupPos[grp_id]
                x_cent = self._fix_pos(grp_pos,x,boxsize=75000.)
                r = AU._dist(x_cent,np.array([0.,0.,0.]))
                r = AU.PhysicalPosition(r,0.19728)

                # r_cut = r <= 150.
                r_cut = r <= 10000.
                rho = rho[r_cut]
                mass = AU.PhysicalMass(np.array(data_dict["MASS"]))[r_cut]
                z = np.array(data_dict["GZ  "])[r_cut]
                u = np.array(data_dict["U   "])[r_cut]
                Nelec = np.array(data_dict["NE  "])[r_cut]
                T = AU.GetTemp(u, Nelec, 5.0/3.0)

                in_ISM = self._find_ISM_gas(rho)
                in_CGM = np.logical_not(in_ISM)
                
                m_ISM = np.append(m_ISM,np.sum(mass[in_ISM]))
                mz_ISM = np.append(mz_ISM,np.sum(z[in_ISM]*mass[in_ISM]))
                m_CGM = np.append(m_CGM,np.sum(mass[in_CGM]))
                mz_CGM = np.append(mz_CGM,np.sum(z[in_CGM]*mass[in_CGM]))
                cool_CGM = np.logical_and(in_CGM,T <= 10.**5.)
                m_cool_CGM = np.append(m_cool_CGM,np.sum(mass[cool_CGM]))
                mz_cool_CGM = np.append(mz_cool_CGM,np.sum(z[cool_CGM]*mass[cool_CGM]))
                warm_CGM = np.logical_and(in_CGM,np.logical_and(T > 10.**5.,T <= 10.**6.))
                m_warm_CGM = np.append(m_warm_CGM,np.sum(mass[warm_CGM]))
                mz_warm_CGM = np.append(mz_warm_CGM,np.sum(z[warm_CGM]*mass[warm_CGM]))
                hot_CGM = np.logical_and(in_CGM,T > 10.**6.)
                m_hot_CGM = np.append(m_hot_CGM,np.sum(mass[hot_CGM]))
                mz_hot_CGM = np.append(mz_hot_CGM,np.sum(z[hot_CGM]*mass[hot_CGM]))

                i += 1

            np.savez(f_npz, m_grp=m_grp,m_winds=m_winds,mz_winds=mz_winds,m_stars=m_stars,mz_stars=mz_stars,\
                m_ISM=m_ISM,mz_ISM=mz_ISM,m_CGM=m_CGM,mz_CGM=mz_CGM,m_cool_CGM=m_cool_CGM,mz_cool_CGM=mz_cool_CGM,\
                m_warm_CGM=m_warm_CGM,mz_warm_CGM=mz_warm_CGM,m_hot_CGM=m_hot_CGM,mz_hot_CGM=mz_hot_CGM,gal_ssfr=gal_ssfr)


        if False:
            # Bin by stellar mass
            n_mbins = 50
            # [mbins_min,mbins_max] = AU._bin_setup(10.**10.,10.**11.5,n_mbins,logbins=True)
            [mbins_min,mbins_max] = AU._bin_setup(10.,11.5,n_mbins)
            mbins_min = 10.**mbins_min
            mbins_max = 10.**mbins_max

            print "mbins_min ",mbins_min

            [m_stars_Q1,m_stars_med,m_stars_Q3] = AU._calc_percentiles_v2(m_stars,mbins_min,mbins_max,m_stars,min_percentile=10,max_percentile=90)
            [m_gal_Q1,m_gal_med,m_gal_Q3] = AU._calc_percentiles_v2(m_stars,mbins_min,mbins_max,m_stars+m_ISM,min_percentile=10,max_percentile=90)
            [mz_gal_Q1,mz_gal_med,mz_gal_Q3] = AU._calc_percentiles_v2(m_stars,mbins_min,mbins_max,mz_stars+mz_ISM,min_percentile=10,max_percentile=90)
            [m_winds_Q1,m_winds_med,m_winds_Q3] = AU._calc_percentiles_v2(m_stars,mbins_min,mbins_max,m_winds,min_percentile=10,max_percentile=90)
            [mz_winds_Q1,mz_winds_med,mz_winds_Q3] = AU._calc_percentiles_v2(m_stars,mbins_min,mbins_max,mz_winds,min_percentile=10,max_percentile=90)
            [m_cool_CGM_Q1,m_cool_CGM_med,m_cool_CGM_Q3] = AU._calc_percentiles_v2(m_stars,mbins_min,mbins_max,m_cool_CGM,min_percentile=10,max_percentile=90)
            [mz_cool_CGM_Q1,mz_cool_CGM_med,mz_cool_CGM_Q3] = AU._calc_percentiles_v2(m_stars,mbins_min,mbins_max,mz_cool_CGM,min_percentile=10,max_percentile=90)
            [m_warm_CGM_Q1,m_warm_CGM_med,m_warm_CGM_Q3] = AU._calc_percentiles_v2(m_stars,mbins_min,mbins_max,m_warm_CGM,min_percentile=10,max_percentile=90)
            [mz_warm_CGM_Q1,mz_warm_CGM_med,mz_warm_CGM_Q3] = AU._calc_percentiles_v2(m_stars,mbins_min,mbins_max,mz_warm_CGM,min_percentile=10,max_percentile=90)
            [m_hot_CGM_Q1,m_hot_CGM_med,m_hot_CGM_Q3] = AU._calc_percentiles_v2(m_stars,mbins_min,mbins_max,m_hot_CGM,min_percentile=10,max_percentile=90)
            [mz_hot_CGM_Q1,mz_hot_CGM_med,mz_hot_CGM_Q3] = AU._calc_percentiles_v2(m_stars,mbins_min,mbins_max,mz_hot_CGM,min_percentile=10,max_percentile=90)
        elif True:
            # Bin by group mass
            n_mbins = 50
            # [mbins_min,mbins_max] = AU._bin_setup(10.**10.,10.**11.5,n_mbins,logbins=True)
            [mbins_min,mbins_max] = AU._bin_setup(11.2,13.,n_mbins)
            mbins_min = 10.**mbins_min
            mbins_max = 10.**mbins_max

            print "mbins_min ",mbins_min

            [m_stars_Q1,m_stars_med,m_stars_Q3] = AU._calc_percentiles_v2(m_grp,mbins_min,mbins_max,m_stars,min_percentile=10,max_percentile=90)
            [m_gal_Q1,m_gal_med,m_gal_Q3] = AU._calc_percentiles_v2(m_grp,mbins_min,mbins_max,m_stars+m_ISM,min_percentile=10,max_percentile=90)
            [mz_gal_Q1,mz_gal_med,mz_gal_Q3] = AU._calc_percentiles_v2(m_grp,mbins_min,mbins_max,mz_stars+mz_ISM,min_percentile=10,max_percentile=90)
            [m_winds_Q1,m_winds_med,m_winds_Q3] = AU._calc_percentiles_v2(m_grp,mbins_min,mbins_max,m_winds,min_percentile=10,max_percentile=90)
            [mz_winds_Q1,mz_winds_med,mz_winds_Q3] = AU._calc_percentiles_v2(m_grp,mbins_min,mbins_max,mz_winds,min_percentile=10,max_percentile=90)
            [m_cool_CGM_Q1,m_cool_CGM_med,m_cool_CGM_Q3] = AU._calc_percentiles_v2(m_grp,mbins_min,mbins_max,m_cool_CGM,min_percentile=10,max_percentile=90)
            [mz_cool_CGM_Q1,mz_cool_CGM_med,mz_cool_CGM_Q3] = AU._calc_percentiles_v2(m_grp,mbins_min,mbins_max,mz_cool_CGM,min_percentile=10,max_percentile=90)
            [m_warm_CGM_Q1,m_warm_CGM_med,m_warm_CGM_Q3] = AU._calc_percentiles_v2(m_grp,mbins_min,mbins_max,m_warm_CGM,min_percentile=10,max_percentile=90)
            [mz_warm_CGM_Q1,mz_warm_CGM_med,mz_warm_CGM_Q3] = AU._calc_percentiles_v2(m_grp,mbins_min,mbins_max,mz_warm_CGM,min_percentile=10,max_percentile=90)
            [m_hot_CGM_Q1,m_hot_CGM_med,m_hot_CGM_Q3] = AU._calc_percentiles_v2(m_grp,mbins_min,mbins_max,m_hot_CGM,min_percentile=10,max_percentile=90)
            [mz_hot_CGM_Q1,mz_hot_CGM_med,mz_hot_CGM_Q3] = AU._calc_percentiles_v2(m_grp,mbins_min,mbins_max,mz_hot_CGM,min_percentile=10,max_percentile=90)
        else:
            # Bin by sSFR:
            self.sub_ids = np.sort(self.sub_ids)
            gal_ssfr = np.array(self.galf['gas_totsfr'])[np.int32(self.sub_ids)]/m_stars
            # impose stellar mass cut 10^11<M<10^11.5
            m_cut = np.logical_and(m_stars >= 10.**11.,m_stars <= 10.**11.5)
            gal_ssfr = gal_ssfr[m_cut]
            m_stars = m_stars[m_cut]
            mz_stars = mz_stars[m_cut]
            m_ISM = m_ISM[m_cut]
            mz_ISM = mz_ISM[m_cut]
            m_cool_CGM = m_cool_CGM[m_cut]
            mz_cool_CGM = mz_cool_CGM[m_cut]
            m_warm_CGM = m_warm_CGM[m_cut]
            mz_warm_CGM = mz_warm_CGM[m_cut]
            m_hot_CGM = m_hot_CGM[m_cut]
            mz_hot_CGM = mz_hot_CGM[m_cut]

            # Bin by ssfr
            n_mbins = 25
            [mbins_min,mbins_max] = AU._bin_setup(-13,-10,n_mbins)
            mbins_min = 10.**mbins_min
            mbins_max = 10.**mbins_max

            print "mbins_min ",mbins_min

            [m_stars_Q1,m_stars_med,m_stars_Q3] = AU._calc_percentiles_v2(gal_ssfr,mbins_min,mbins_max,m_stars,min_percentile=10,max_percentile=90)
            [m_gal_Q1,m_gal_med,m_gal_Q3] = AU._calc_percentiles_v2(gal_ssfr,mbins_min,mbins_max,m_stars+m_ISM,min_percentile=10,max_percentile=90)
            [mz_gal_Q1,mz_gal_med,mz_gal_Q3] = AU._calc_percentiles_v2(gal_ssfr,mbins_min,mbins_max,mz_stars+mz_ISM,min_percentile=10,max_percentile=90)
            [m_winds_Q1,m_winds_med,m_winds_Q3] = AU._calc_percentiles_v2(gal_ssfr,mbins_min,mbins_max,m_winds,min_percentile=10,max_percentile=90)
            [mz_winds_Q1,mz_winds_med,mz_winds_Q3] = AU._calc_percentiles_v2(gal_ssfr,mbins_min,mbins_max,mz_winds,min_percentile=10,max_percentile=90)
            [m_cool_CGM_Q1,m_cool_CGM_med,m_cool_CGM_Q3] = AU._calc_percentiles_v2(gal_ssfr,mbins_min,mbins_max,m_cool_CGM,min_percentile=10,max_percentile=90)
            [mz_cool_CGM_Q1,mz_cool_CGM_med,mz_cool_CGM_Q3] = AU._calc_percentiles_v2(gal_ssfr,mbins_min,mbins_max,mz_cool_CGM,min_percentile=10,max_percentile=90)
            [m_warm_CGM_Q1,m_warm_CGM_med,m_warm_CGM_Q3] = AU._calc_percentiles_v2(gal_ssfr,mbins_min,mbins_max,m_warm_CGM,min_percentile=10,max_percentile=90)
            [mz_warm_CGM_Q1,mz_warm_CGM_med,mz_warm_CGM_Q3] = AU._calc_percentiles_v2(gal_ssfr,mbins_min,mbins_max,mz_warm_CGM,min_percentile=10,max_percentile=90)
            [m_hot_CGM_Q1,m_hot_CGM_med,m_hot_CGM_Q3] = AU._calc_percentiles_v2(gal_ssfr,mbins_min,mbins_max,m_hot_CGM,min_percentile=10,max_percentile=90)
            [mz_hot_CGM_Q1,mz_hot_CGM_med,mz_hot_CGM_Q3] = AU._calc_percentiles_v2(gal_ssfr,mbins_min,mbins_max,mz_hot_CGM,min_percentile=10,max_percentile=90)

        if BH_split:
            # In each group mass bin, take separately the groups with top 25% BH mass, and groups with bot 25% BH mass [NOT IMPLEMENTED]
            # Do simple thing: split by arbitrary BH mass:
            gal_BHmass = AU.PhysicalMass(np.array(self.galf['bh_totmass'])[np.int32(self.sub_ids)])





        # Plot Baryon budget
        plt.figure(figsize=(5,5))
        x = np.log10(mbins_min) 
        print "x ",x
        # plt.plot(x,np.log10(m_winds_med),label='Winds',color='brown',ls='dotted')
        # plt.fill_between(x,np.log10(m_winds_Q1),np.log10(m_winds_Q3),alpha=0.3,color='brown')
        plt.plot(x,np.log10(m_gal_med),label='Stars+ISM',color='green')
        # plt.fill_between(x,np.log10(m_gal_Q1),np.log10(m_gal_Q3),alpha=0.3,color='green')
        plt.plot(x,np.log10(m_cool_CGM_med),label='Cool CGM',color='blue')
        # plt.fill_between(x,np.log10(m_cool_CGM_Q1),np.log10(m_cool_CGM_Q3),alpha=0.3,color='blue')
        plt.plot(x,np.log10(m_warm_CGM_med),label='Warm CGM',color='gold')
        # plt.fill_between(x,np.log10(m_warm_CGM_Q1),np.log10(m_warm_CGM_Q3),alpha=0.3,color='gold')
        plt.plot(x,np.log10(m_hot_CGM_med),label='Hot CGM',color='red')
        # plt.fill_between(x,np.log10(m_hot_CGM_Q1),np.log10(m_hot_CGM_Q3),alpha=0.3,color='red')
        plt.xlabel("Log10(Halo Mass)")
        plt.ylabel("Log10 (Mass)")
        plt.legend(loc=4)


        # x = np.ones(5)*12.2
        x = np.array([12.1,12.15,12.2,12.25,12.3])
        y = np.array([10.7,11.0,10.5,9.6,11.3])
        yerr = np.array([0.1,0.08,0.5,0.6,0.2])
        plt.errorbar(np.copy(x[1]),np.copy(y[0]),yerr=np.copy(yerr[0]),color='green',lw=6)
        plt.errorbar(np.copy(x[0]),np.copy(y[1]),yerr=np.copy(yerr[1]),color='blue',lw=6)
        plt.errorbar(np.copy(x[2]),np.copy(y[2]),yerr=np.copy(yerr[2]),color='gold',lw=6)
        plt.errorbar(np.copy(x[3]),np.copy(y[3]),yerr=np.copy(yerr[3]),color='red',lw=6)
        # x = np.array([10.,15.])
        # plt.fill_between(x,np.log10(4*10**10.),np.log10(7*10**10.),color='green',alpha=0.2)
        # plt.fill_between(x,np.log10(7*10**10.),np.log10(12*10**10.),color='blue',alpha=0.2)
        # plt.fill_between(x,np.log10(1*10**10.),np.log10(10*10**10.),color='gold',alpha=0.2)
        # plt.fill_between(x,np.log10(1*10**9.),np.log10(14*10**9.),color='red',alpha=0.2)
        # plt.errorbar(np.copy(x[4]),np.copy(y[4]),yerr=np.copy(yerr[4]),color='black')
        # plt.xlim([11.,13.])
        plt.savefig(self.fig_base+'mass_budget_grp.pdf', bbox_inches='tight') 
        

        # Plot metal budget
        plt.figure() 
        x = np.log10(mbins_min)
        # plt.plot(x,np.log10(mz_winds_med),label='Winds',color='brown',ls='dotted')
        # plt.fill_between(x,np.log10(mz_winds_Q1),np.log10(mz_winds_Q3),alpha=0.3,color='brown')
        plt.plot(x,np.log10(mz_gal_med),label='Stars+ISM',color='green')
        # plt.fill_between(x,np.log10(mz_gal_Q1),np.log10(mz_gal_Q3),alpha=0.3,color='green')
        plt.plot(x,np.log10(mz_cool_CGM_med),label='Cool CGM',color='blue')
        # plt.fill_between(x,np.log10(mz_cool_CGM_Q1),np.log10(mz_cool_CGM_Q3),alpha=0.3,color='blue')
        plt.plot(x,np.log10(mz_warm_CGM_med),label='Warm CGM',color='gold')
        # plt.fill_between(x,np.log10(mz_warm_CGM_Q1),np.log10(mz_warm_CGM_Q3),alpha=0.3,color='gold')
        plt.plot(x,np.log10(mz_hot_CGM_med),label='Hot CGM',color='red')
        # plt.fill_between(x,np.log10(mz_hot_CGM_Q1),np.log10(mz_hot_CGM_Q3),alpha=0.3,color='red')
        plt.plot(x,np.log10(mz_cool_CGM_med+mz_warm_CGM_med+mz_hot_CGM_med),label='CGM',color='black')
        plt.savefig(self.fig_base+'metal_budget_grp.pdf')





    def OVI_weighted_gas(self,data_type,Mmin=10**11.4,Mmax=10**11.45,low_ssfr_pop=False,high_ssfr_pop=False,savename=None):
        # Plot properties of OVI-weighted gas, vs all gas.

        plt.close('all')
        plt.figure(figsize=(3.54331,3.14))

        # Set up bins for radius:
        n_Rbins = 50
        [Rbins_min,Rbins_max] = AU._bin_setup(0.,300.,n_Rbins)
        Rbins_med = (Rbins_min+Rbins_max)/2.

        # pseudo code:
        # find all groups in central galaxy stellar mass/sSFR bin that we want

        # get desired group ids:
        # grp_ids = np.arange(np.size(self.gal_sm))
        # mass_cut = np.logical_and(self.gal_sm > Mmin, self.gal_sm < Mmax)
        # grp_ids = grp_ids[mass_cut]
        # if low_ssfr_pop:
        #     save_ssfr = self.gal_ssfr[mass_cut]
        #     ssfr_cut = save_ssfr < 10.**-11
        #     grp_ids = grp_ids[ssfr_cut]
        # elif high_ssfr_pop:
        #     save_ssfr = self.gal_ssfr[mass_cut]
        #     ssfr_cut = save_ssfr >= 10.**-11
        #     grp_ids = grp_ids[ssfr_cut]
        # selecting by myself:
        grp_ids = np.array([19])
        n_grps = np.size(grp_ids)

        print "grp_ids ",grp_ids

        z_profile_gal = np.zeros([n_grps,n_Rbins])
        T_profile_gal = np.zeros([n_grps,n_Rbins])
        rho_profile_gal = np.zeros([n_grps,n_Rbins])
        z_OVI_gal = np.zeros([n_grps,n_Rbins])
        T_OVI_gal = np.zeros([n_grps,n_Rbins])
        rho_OVI_gal = np.zeros([n_grps,n_Rbins])

        # load CGM snapshots for these groups
        for j in np.arange(np.size(grp_ids)):
            grp_id = grp_ids[j]
            fn = self.CGMsnap_base + "s{}/{}.hdf5".format(self.snapnum,str(int(grp_id)).zfill(5))
            print "fn ",fn
            grp_data = self.load_CGM_snap(fn)
            m = AU.PhysicalMass(np.array(grp_data['MASS']))
            z = np.array(grp_data["GZ  "])
            met = np.array(grp_data["GMET"])
            T = self._get_T(grp_data)
            rho = AU.PhysicalDensity(np.array(grp_data["RHO "]),self.redshift)
            r = AU.PhysicalPosition(self._calc_radii(grp_data),self.redshift)

            # Get OVI fraction using CLOUDY table
            rho_Hatoms = rho*(met[:,0]/AU.ProtonMass)
            # Cloudy cut:
            rho_cut = np.logical_and(np.log10(rho_Hatoms) >= -7.,np.log10(rho_Hatoms) <= 4.)
            T_cut = np.logical_and(np.log10(T) >= 3.,np.log10(T) <= 8.6)
            cloudy_cut = np.logical_and(T_cut,rho_cut)
            if np.sum(cloudy_cut) < np.size(m): 
                print "cloudy cut removed {} particles".format(np.size(m)-np.sum(cloudy_cut))
                r = r[cloudy_cut]
                m = m[cloudy_cut]
                met = met[cloudy_cut]
                z = z[cloudy_cut]
                T = T[cloudy_cut]
                rho = rho[cloudy_cut]
                rho_Hatoms = rho_Hatoms[cloudy_cut]
            beta = self.gal_ssfr[grp_id]/(r**2.)
            elem_frac = met[:,AU.elem_lookup('O')]
            ion_frac = self.tab.ion('O',6,rho_Hatoms,T,beta=beta)
            OVI_frac = elem_frac*ion_frac
            OVI_mass = m*OVI_frac

            
            
            # Calculate mass-weighted rho, T as a function of radius:

            for k in np.arange(n_Rbins):
                in_Rbin = np.logical_and(r > Rbins_min[k], r < Rbins_max[k])
                z_profile_gal[j,k] = np.sum(m[in_Rbin]*z[in_Rbin])/np.sum(m[in_Rbin])
                T_profile_gal[j,k] = np.sum(m[in_Rbin]*T[in_Rbin])/np.sum(m[in_Rbin])
                rho_profile_gal[j,k] = np.sum(m[in_Rbin]*rho[in_Rbin])/np.sum(m[in_Rbin])
                # Also calculate OVI-weighted values...
                z_OVI_gal[j,k] = np.sum(OVI_mass[in_Rbin]*z[in_Rbin])/np.sum(OVI_mass[in_Rbin])
                T_OVI_gal[j,k] = np.sum(OVI_mass[in_Rbin]*T[in_Rbin])/np.sum(OVI_mass[in_Rbin])
                rho_OVI_gal[j,k] = np.sum(OVI_mass[in_Rbin]*rho[in_Rbin])/np.sum(OVI_mass[in_Rbin])


        if data_type == 'z':
            [Q1,med,Q3] = AU._calc_percentiles(z_profile_gal)
            [Q1,med_OVI,Q3] = AU._calc_percentiles(z_OVI_gal)
            plt.ylabel(r'$\log_{10}\left[\frac{Z}{Z_\odot}\right]$')
            med /= AU.SolarMetallicity
            med_OVI /= AU.SolarMetallicity
        if data_type == 'T':
            [Q1,med,Q3] = AU._calc_percentiles(T_profile_gal)
            [Q1,med_OVI,Q3] = AU._calc_percentiles(T_OVI_gal)
            plt.ylabel(r'$\log_{10}\left[\frac{{\rm Temperature}}{\rm K}\right]$')
        if data_type == 'rho':
            [Q1,med,Q3] = AU._calc_percentiles(rho_profile_gal)
            [Q1,med_OVI,Q3] = AU._calc_percentiles(rho_OVI_gal)
            plt.ylabel(r'Density [cm$^{-3}$]')

        plt.plot(Rbins_min,np.log10(med),label='All Gas')
        plt.plot(Rbins_min,np.log10(med_OVI),label='OVI Gas')
        plt.xlabel(r'3D Radius [R$_{200}$]')
        plt.legend()

        if savename != None:
            plt.savefig(self.fig_base+savename+".pdf")
        else:
            plt.savefig(self.fig_base+"{}_radprof.pdf".format(data_type))




    def galprop_vs_CGMprop(self,galprop,CGMprop,savename=None):
        # galprops: stellar mass, SFR, environment
        # CGMprop: total metallicity, total metal mass, ion covering fraction

        # Get CGMsnap group ids, and main subhalo ids:
        self.load_CGMsnap_ids()
        cat = readsubfHDF5.subfind_catalog('/n/ghernquist/Illustris/Runs/Illustris-1/',self.snapnum,keysel=['GroupFirstSub'],subcat=False)
        self.sub_ids = cat.GroupFirstSub[np.int32(self.grp_ids)]

        print "Loading galaxy properties... "
        if galprop == 'sm':
            galf = h5py.File(self.snapbase+"/postprocessing/galprop/galprop_{}.hdf5".format(self.snapnum),'r')
            gal_sm = AU.PhysicalMass(np.array(galf['stellar_totmass']))
            x = gal_sm[self.sub_ids]
        elif galprop == 'SFR':
            galf = h5py.File(self.snapbase+"/postprocessing/galprop/galprop_{}.hdf5".format(self.snapnum),'r')
            gal_sfr = AU.PhysicalMass(np.array(galf['gas_totsfr']))
            x = gal_sfr[self.sub_ids]
        elif galprop == 'sSFR':
            galf = h5py.File(self.snapbase+"/postprocessing/galprop/galprop_{}.hdf5".format(self.snapnum),'r')
            gal_sm = AU.PhysicalMass(np.array(galf['stellar_totmass']))
            gal_sfr = AU.PhysicalMass(np.array(galf['gas_totsfr']))
            x = gal_sfr[self.sub_ids]/gal_sm[self.sub_ids]


        print "Done loading galaxy properties! "

        print "Loading CGM properties... "
        # def load_CGM_snap(fn,load='all'):
        #     data_dict = {}
        #     f = h5py.File(fn,'r')
        #     if load=='all':
        #         for key in f['Header'].attrs:
        #             data_dict[key] = f['Header'].attrs[key]
        #         for key in f['PartType0']:
        #             data_dict[key] = np.copy(f['PartType0'][key])
        #     else:
        #         data_dict[key] = f['Header'].attrs[key]
        #     return data_dict


        if CGMprop == 'z_with_ISM':
            # Total metal mass including ISM

            # First check if it has been written to an npz already:
            f_npz = self.npz_base + "{}.npz".format(CGMprop)
            if os.path.isfile(f_npz):
                dat = np.load(f_npz)
                zmass = f_npz['zmass']
            else:
                print "File {} does not exist!  Will write it after loading...".format(f_npz)
            zmass = np.zeros(0)
            i = 0
            for fn in glob.glob(self.CGMsnap_base+"s{}/*.hdf5".format(self.snapnum)):
                print "i ",i
                data_dict = load_CGM_snap(fn)
                zmass = np.append(zmass,AU.PhysicalMass(np.sum(np.array(data_dict["GZ  "])*np.array(data_dict["MASS"]))))
                i += 1
            y = zmass
            np.savez(f_npz,zmass=zmass)


        elif CGMprop == 'z_without_ISM': 
            # First check if it has been written to an npz already:
            f_npz = self.npz_base + "{}.npz".format(CGMprop)
            if os.path.isfile(f_npz):
                dat = np.load(f_npz)
                CGM_mass = dat['CGM_mass']
                ISM_mass = dat['ISM_mass']
                CGM_zmass = dat['CGM_zmass']
                ISM_zmass = dat['ISM_zmass']
            else:
                print "File {} does not exist!  Will write it after loading...".format(f_npz)
                CGM_mass = np.zeros(0)
                ISM_mass = np.zeros(0)
                CGM_zmass = np.zeros(0)
                ISM_zmass = np.zeros(0)
                i = 0
                for fn in glob.glob(self.CGMsnap_base+"s{}/*.hdf5".format(self.snapnum)):
                    print "i ",i
                    data_dict = load_CGM_snap(fn)
                    rho = AU.PhysicalDensity(np.array(data_dict["RHO "]),0.19728)
                    in_ISM = self._find_ISM_gas(rho)
                    in_CGM = np.logical_not(in_ISM)
                    mass = AU.PhysicalMass(np.array(data_dict["MASS"]))
                    z = np.array(data_dict["GZ  "])
                    CGM_mass = np.append(CGM_mass,np.sum(mass[in_CGM]))
                    ISM_mass = np.append(ISM_mass,np.sum(mass[in_ISM]))
                    CGM_zmass = np.append(CGM_zmass,np.sum(z[in_CGM]*mass[in_CGM]))
                    ISM_zmass = np.append(ISM_zmass,np.sum(z[in_ISM]*mass[in_ISM]))
                    i += 1
                np.savez(f_npz,CGM_mass=CGM_mass,ISM_mass=ISM_mass,CGM_zmass=CGM_zmass,ISM_zmass=ISM_zmass)
        
        elif CGMprop == 'metallicity':
             # First check if it has been written to an npz already:
            f_npz = self.npz_base + "{}.npz".format("z_without_ISM")
            if os.path.isfile(f_npz):
                dat = np.load(f_npz)
                CGM_mass = dat['CGM_mass']
                ISM_mass = dat['ISM_mass']
                CGM_zmass = dat['CGM_zmass']
                ISM_zmass = dat['ISM_zmass']
                z_CGM = CGM_zmass/CGM_mass
                z_ISM = ISM_zmass/ISM_mass
            else:
                print "File {} does not exist!  Halp...".format(f_npz)

        elif CGMprop == 'CGM_ISM_metal_ratio':
             # First check if it has been written to an npz already:
            f_npz = self.npz_base + "{}.npz".format("z_without_ISM")
            if os.path.isfile(f_npz):
                dat = np.load(f_npz)
                CGM_mass = dat['CGM_mass']
                ISM_mass = dat['ISM_mass']
                CGM_zmass = dat['CGM_zmass']
                ISM_zmass = dat['ISM_zmass']
            else:
                print "File {} does not exist!  Halp...".format(f_npz)


        print "Done loading CGM properties!"

        plt.figure()
        # plt.scatter(np.log10(x),np.log10(CGM_mass),marker='.',s=10,alpha=0.3,label='CGM_mass')
        # plt.scatter(np.log10(x),np.log10(ISM_mass),marker='.',s=10,alpha=0.3,label='ISM_mass',color='green')
        # plt.scatter(np.log10(x),np.log10(CGM_zmass),marker='.',s=10,alpha=0.3,label='CGM_zmass',color='red')
        # plt.scatter(np.log10(x),np.log10(ISM_zmass),marker='.',s=10,alpha=0.3,label='ISM_zmass',color='magenta')
        # plt.ylim([7.,12.])
        # plt.scatter(np.log10(x),np.log10(z_CGM),marker='.',s=10,alpha=0.3,label='z_CGM')
        # plt.scatter(np.log10(x),np.log10(z_ISM),marker='.',s=10,alpha=0.3,label='z_ISM',color='green')
        plt.scatter(np.log10(x),np.log10(CGM_zmass/ISM_zmass),marker='.',s=10,alpha=0.3,label='metal ratio')
        # plt.legend()
        plt.savefig(self.fig_base+"mass_vs_metalratio.pdf", bbox_inches='tight')

        # if savename != None: 
        #     plt.savefig(self.fig_base+savename+".pdf", bbox_inches='tight')
        # else:
        #     plt.savefig(self.fig_base+"galprop_CGMprop_test.pdf", bbox_inches='tight')

        

        def load_all_CGMsnaps():
            data_dict = {}
            for fn in glob.glob(self.CGMsnap_base+"s{}/*.hdf5".format(self.snapnum)):
                print "fn ",fn
                f = h5py.File(fn,'r')
                grp_id = f['Header'].attrs['grp_id']
                for key in f['grids']:
                    data_dict[key] = f['grids'][key]
                grid_rad = f['Header'].attrs['grid_radius_pkpc']
                ngrid = f['Header'].attrs['ngrid']
                grid = np.array(f['grids'][species])
                f.close()

                self.snapnum



    def coldens_plot(self,species,kpc_mode=False,Rvir_mode=False,vmin=13.,vmax=16.,low_ssfr_pop=False,high_ssfr_pop=False,Mmin=10**10.5,Mmax=10**11.,cloudy=None,tempfac=None,savename=None):
        # Load all grid files.
        # Calculate radii for each one.
        # hexbin the full radius vs column density spread, by kpc or Rvir.

        nbinx = 100.
        nbiny = 100.
        twod_hist = np.zeros([nbinx,nbiny])


        
        # save as hdf5 file.
        if Mmin == 10.**11:
            masstag = 'highm'
        elif Mmin == 10.**10.5:
            masstag = 'medm'
        elif Mmin == 10.**10.:
            masstag = 'lowm'
        if low_ssfr_pop: ssfr_tag = 'lowssfr'
        elif high_ssfr_pop: ssfr_tag = 'highssfr'
        if tempfac != None and cloudy != None:
            save_hdf5 = self.npz_base+"{}_{}_{}_{}_t{}.hdf5".format(species,masstag,ssfr_tag,cloudy,tempfac)
        else: 
            save_hdf5 = self.npz_base+"{}_{}_{}_{}.hdf5".format(species,masstag,ssfr_tag,cloudy)
        
        print "save_hdf5 ",save_hdf5

        if os.path.isfile(save_hdf5):
            g = h5py.File(save_hdf5,'r')
            twod_hist = np.array(g['histogram']['twod_hist'])
            xedges = np.array(g['histogram']['xedges'])
            yedges = np.array(g['histogram']['yedges'])
            g.close()
        else: 
            i = 0
            print "hi"
            print self.grid_base+"s{}/*.hdf5".format(self.snapnum)
            print glob.glob(self.grid_base+"s{}/*.hdf5".format(self.snapnum))
            for fn in glob.glob(self.grid_base+"s{}/*.hdf5".format(self.snapnum)):
                print "fn ",fn
                print "i ",i
                f = h5py.File(fn,'r')
                grp_id = f['Header'].attrs['grp_id']
                grp_Rvir = f['Header'].attrs['grp_Rvir']
                grid_rad = f['Header'].attrs['grid_radius_pkpc']
                ngrid = f['Header'].attrs['ngrid']
                if tempfac != None and cloudy != None:
                    grid = np.array(f[cloudy]['temp_fac_{}'.format(tempfac)][species])
                elif tempfac == None and cloudy != None:
                    grid = np.array(f[cloudy][species])
                else:
                    grid = np.array(f['grids'][species])
                f.close()

                # save as hdf5 file?

                # check if galaxy matches desired properties.
                gal_props = self._get_gal_props(grp_id)
                ssfr = gal_props['ssfr']
                sm = gal_props['sm']

                cond1 = (ssfr < 10.**-11.) and (sm > Mmin) and (sm < Mmax) #(sm > 10.**10.5) and (sm < 10.**11)
                cond2 = (ssfr > 10.**-11.) and (sm > Mmin) and (sm < Mmax) #(sm > 10.**10.5) and (sm < 10.**11)
                if (low_ssfr_pop and cond1) or (high_ssfr_pop and cond2) or (not low_ssfr_pop and not high_ssfr_pop):
                    [gridx,gridy] = np.meshgrid(np.arange(ngrid),np.arange(ngrid))
                    grid_cent = (ngrid-1)/2. #assume square grid: grid_centx = grid_centy = grid_cent
                    r_grid = np.sqrt((gridx-grid_cent)**2+(gridy-grid_cent)**2)
                    r_kpc = self._grid_to_kpc(r_grid,ngrid,grid_rad)
                    r_Rvir = r_kpc/AU.PhysicalPosition(grp_Rvir,0.19728)

                    if Rvir_mode:
                        H,xedges,yedges = np.histogram2d(np.ravel(grid),np.ravel(r_Rvir),bins=[nbinx,nbiny],range=[[0.,1.],[vmin,vmax]])
                        twod_hist += H
                    elif kpc_mode:
                        # H,xedges,yedges = np.histogram2d(np.ravel(r_kpc),np.ravel(grid),bins=[nbinx,nbiny],range=[[0.,150.],[vmin,vmax]])
                        H,yedges,xedges = np.histogram2d(np.ravel(grid),np.ravel(r_kpc),bins=[nbiny,nbinx],range=[[vmin,vmax],[0.,200.]],normed=False)
                        twod_hist += H
                    i+=1

            # g = h5py.File(save_hdf5,'w-')
            # h_grp = g.create_group('histogram')
            # h_grp.create_dataset("twod_hist",data=twod_hist)
            # h_grp.create_dataset("xedges",data=xedges)
            # h_grp.create_dataset("yedges",data=yedges)
            # g.close()

        # np.savez(f_npz,N=N,r_kpc=r_kpc,r_Rvir=r_Rvir)




        # a = np.array([[1,2],[3,4]])
        # plt.imshow(a,aspect='auto',cmap=plt.cm.cubehelix,origin="lower",extent=(0,150.,vmin,vmax))
        print "twod_hist ",twod_hist
        print "xedges ",xedges
        for i in np.arange(nbinx):
            twod_hist[i] /= (2*np.pi*(xedges[1:]))
        plt.figure(figsize=(4,4))
        plt.imshow(twod_hist,aspect='auto',cmap=plt.cm.cubehelix,origin="lower",extent=(0,200.,vmin,vmax),zorder=1)
        plt.ylabel(r"log$_{10}$ N$_\mathrm{OVI}$ (cm$^{-2}$)")
        plt.xlabel("Projected Radius (pkpc)")
        
        # Add COS-Halos stuff:
        cos_dat = np.loadtxt("cosdata_condense.txt")
        elem = cos_dat[:,0]
        ion = cos_dat[:,1]
        M = 10.**(cos_dat[:,2])
        sfr = cos_dat[:,3]
        sfr_upperlim = cos_dat[:,4]
        R = cos_dat[:,5]
        N_spec = cos_dat[:,6]
        N_spec_err = cos_dat[:,7]
        
        ssfr = sfr/M
        if low_ssfr_pop:
            ssfr_cut = ssfr <= 10.**-11.
            if species == "Mg2": species_cut = np.logical_and(elem==6,ion==2)
            elif species == "Si3": species_cut = np.logical_and(elem==7,ion==3)
            elif species == "O6": species_cut = np.logical_and(elem==4,ion==6)
            elif species == "N5": species_cut = np.logical_and(elem==3,ion==5)

            mass_cut = np.logical_and(M>=Mmin,M<=Mmax) #M>=10.**11. 
            full = np.logical_and(np.logical_and(ssfr_cut,species_cut),mass_cut)
            x = R[full]
            y = N_spec[full]
            yerr = N_spec_err[full]
            upper_lim = yerr == -2
            normal = np.logical_not(upper_lim)
            print " x[normal],y[normal],yerr[normal]",x[normal],y[normal],yerr[normal]
            plt.errorbar(x[normal],y[normal],yerr=yerr[normal],ls='none',marker='o',color='red',zorder=2)
            plt.errorbar(x[upper_lim],y[upper_lim],ls='none',marker='v',color='red',zorder=3)
        elif high_ssfr_pop:
            ssfr_cut = ssfr >= 10.**-11.
            if species == "Mg2": species_cut = np.logical_and(elem==6,ion==2)
            elif species == "Si3": species_cut = np.logical_and(elem==7,ion==3)
            elif species == "O6": species_cut = np.logical_and(elem==4,ion==6)
            mass_cut = np.logical_and(M>=Mmin,M<=Mmax) #M>=10.**11. 
            full = np.logical_and(np.logical_and(ssfr_cut,species_cut),mass_cut)
            x = R[full]
            y = N_spec[full]
            yerr = N_spec_err[full]
            upper_lim = yerr == -2
            normal = np.logical_not(upper_lim)
            print " x[normal],y[normal],yerr[normal]",x[normal],y[normal],yerr[normal]
            plt.errorbar(x[normal],y[normal],yerr=yerr[normal],ls='none',marker='o',color='cyan',zorder=2)
            plt.errorbar(x[upper_lim],y[upper_lim],ls='none',marker='v',color='cyan',zorder=3)


        if species == 'N5':
            x_ul = np.array([78.,82.9,76.9,112.2,134.6,101.0,83.1,91.0,110.3,120.9,149.2,90.7,60.3,43.7,23.0,52.7,18.3,54.7,19.1,154.4,132.3,31.6,112.1,37.2,87.2,35.4,32.3,88.3,82.7,38.9,142.8,113.1,47.1,102.5,150.0,33.8,116.4,46.2])
            N_ul = np.array([13.5,13.6,13.5,13.4,13.6,13.4,13.7,14.1,13.9,13.5,13.5,13.8,14.5,14.0,14.0,13.8,13.8,13.5,13.6,13.7,13.4,13.8,13.7,14.0,13.5,13.9,13.5,13.7,14.0,14.0,14.0,13.7,13.7,13.7,13.7,13.9,13.6,13.5])
            x_det = np.array([44.3,19.7,92.6,96.8])
            N_det = np.array([14.2,13.7,13.7,13.7])
            plt.errorbar(x_ul,N_ul,ls='none',marker='v',color='yellow')
            plt.errorbar(x_det,N_net,yerr=np.ones_like(x_det)*0.1,ls='none',marker='o',color='green')

            # h = np.loadtxt('/n/home04/jsuresh/Python/MyModules/CGM/coshalos_NV.txt')
            # h =

        plt.ylim([vmin,vmax])
        if savename == None:
            plt.savefig(self.fig_base+"{}_coldens_medm_lowssfr_wCOS".format(species), bbox_inches='tight')
        else:
            plt.savefig(self.fig_base+savename+".pdf",bbox_inches='tight')



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


        for fn in glob.glob(self.grid_base+"s{}/*.hdf5".format(self.snapnum)):
            print "fn ",fn
            f = h5py.File(fn,'r')
            grp_id = f['Header'].attrs['grp_id']
            grid_rad = f['Header'].attrs['grid_radius_pkpc']
            ngrid = f['Header'].attrs['ngrid']
            # grid = np.array(f['grids'][species])
            # for cloudy in ['ion_out_fancy_atten','UVB_sf_xrays_ext']:
            for cloudy in ['UVB_sf_xrays_ext']:
                grid = np.array(f[cloudy][species])

                img_savepath = "{}/grids/s{}/{}_{}_{}_mult.pdf".format(self.fig_base,self.snapnum,str(int(grp_id)).zfill(5),species,cloudy)

                maxdist = grid_rad

                plt.close('all')
                if HI_custom_map:
                    plt.imshow(grid,origin='lower',extent=(-maxdist,maxdist,-maxdist,maxdist),vmin=vmin,vmax=vmax,cmap=my_cmap) 
                else: 
                    plt.imshow(grid,origin='lower',extent=(-maxdist,maxdist,-maxdist,maxdist),vmin=vmin,vmax=vmax,cmap=plt.cm.cubehelix) 
                bar=plt.colorbar()
                bar_label = r"log$_{10}$ N$_\mathrm{"+species+"}$ (cm$^{-2}$)"
                bar.set_label(bar_label)
                plt.xlabel(r"y (pkpc)")
                plt.ylabel(r"z (pkpc)")
                plt.savefig(img_savepath)
            f.close()


    def load_CGM_snap(self,fn,load='all'):
        data_dict = {}
        # Check that file exists:
        if os.path.isfile(fn):
            f = h5py.File(fn,'r')
            if load=='all':
                for key in f['Header'].attrs:
                    data_dict[key] = f['Header'].attrs[key]
                for key in f['PartType0']:
                    data_dict[key] = np.copy(f['PartType0'][key])
            else:
                data_dict[key] = f['Header'].attrs[key]
            return data_dict
        else:
            print "fn doesnt exist "

    def load_CGM_snap_ids(self):
        grp_ids = np.zeros(0)
        i = 0
        for fn in glob.glob(self.CGMsnap_base+"s{}/*.hdf5".format(self.snapnum)):
            f = h5py.File(fn,'r')
            grp_ids = np.append(grp_ids,f['Header'].attrs['grp_id'])
            f.close()
            i += 1
        self.grp_ids = np.sort(grp_ids)

    # def find_all_grid_files(self):
    # def load_grid_data(self):


    def load_gal_props(self):
        galf = h5py.File(self.snapbase+"/postprocessing/galprop/galprop_{}.hdf5".format(self.snapnum),'r')
        gal_sm = np.array(galf['stellar_totmass'])
        gal_sfr = np.array(galf['gas_totsfr'])

        sm_mass_select = np.logical_and(gal_sm > 1.,gal_sm < 40)
        gal_ids = np.arange(np.size(gal_sm))[sm_mass_select]

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
        grnr = self.cat.SubhaloGrNr[primary_gal_ids]
        if np.size(np.unique(grnr)) < np.size(grnr): print "had non-uniques!"

        self.gal_sm = AU.PhysicalMass(gal_sm[primary_gal_ids])
        self.gal_sfr = gal_sfr[primary_gal_ids]
        self.gal_ssfr = self.gal_sfr/self.gal_sm


    def _find_ISM_gas(self,rho):
        m_typ = 2.1e-24 # mass of typical particle (75% H, 25% He)
        SF = rho/m_typ > 0.13
        return SF

    def _grid_to_kpc(self,grid_dist,ngrid,grid_radius):
        kpc_dist = grid_dist * (2*grid_radius)/ngrid
        return kpc_dist

    def _get_gal_props(self,grp_id):
        sub_id = np.int32(self.cat.GroupFirstSub[grp_id])
        
        data_dict = {}
        data_dict['sm'] = AU.PhysicalMass(np.array(self.galf['stellar_totmass'][sub_id]))
        data_dict['sfr'] = np.array(self.galf['gas_totsfr'][sub_id])
        data_dict['ssfr'] = data_dict['sfr']/data_dict['sm']
        return data_dict


    def _fix_pos(self,grp_pos,pos,boxsize):

        def _pbc(delta):
            index1=delta > +boxsize/2
            index2=delta < -boxsize/2
            delta[index1]-=boxsize
            delta[index2]+=boxsize
            return delta

        dx = _pbc(pos[:,0]-grp_pos[0])
        dy = _pbc(pos[:,1]-grp_pos[1])
        dz = _pbc(pos[:,2]-grp_pos[2])

        new_pos = np.zeros_like(pos)
        new_pos[:,0] = dx
        new_pos[:,1] = dz
        new_pos[:,2] = dy

        return new_pos


    def load_CGM_snap(self,fn,load='all'):
        data_dict = {}
        f = h5py.File(fn,'r')
        if load=='all':
            for key in f['Header'].attrs:
                data_dict[key] = f['Header'].attrs[key]
            for key in f['PartType0']:
                data_dict[key] = np.copy(f['PartType0'][key])
        else:
            data_dict[key] = f['Header'].attrs[key]
        return data_dict

    def _calc_radii(self,grp_data):
        # box = grp_data['box']
        print "Inputting box directly"
        box = 75000.
        grp_pos = grp_data['grp_pos']
        pos = grp_data['POS ']
        r = AU._pbc_dist(pos,grp_pos,boxsize=box)
        return r

    def _get_T(self,grp_data):
        u = np.array(grp_data['U   '])#'InternalEnergy'
        nelec = np.array(grp_data['NE  '])# "ElectronAbundance"
        T = AU.GetTemp(u, nelec, gamma = 5.0/3.0)
        return T

if __name__ == '__main__':
    illustris_fan()
