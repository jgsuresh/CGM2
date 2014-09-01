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

import glob
import os

class illustris_fan:
    def __init__(self):
        self.snapbase = "/n/ghernquist/Illustris/Runs/Illustris-1/output/"
        self.CGMsnap_base = '/n/home04/jsuresh/scratch1/AREPOfest/data/CGM_snaps/'
        self.grid_base = '/n/home04/jsuresh/scratch1/AREPOfest/data/grids/'
        self.fig_base = '/n/home04/jsuresh/scratch1/AREPOfest/data/figs/'
        self.npz_base = '/n/home04/jsuresh/scratch1/AREPOfest/data/npz/'

        self.snapnum = 120

        #read in catalog
        # self.cat = readsubfHDF5.subfind_catalog(self.snapbase, self.snapnum, long_ids=True, double_output=True, keysel=["GroupFirstSub","SubhaloGrNr"])
        # self.load_gal_props()
        # self.gal_mass_vs_sSFR()

        # self.plot_grids("H1",vmax=25.)
        # self.plot_grids("SiIII",vmax=16.)
        # self.plot_grids("O6",vmax=16.)
        # self.galprop_vs_CGMprop('sm','CGM_ISM_metal_ratio')
        # self.gal_mass_vs_Rvir()
        self.coldens_plot('O6',kpc_mode=True)


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
        cat = readsubfHDF5.subfind_catalog('/n/ghernquist/Illustris/Runs/Illustris-1/',self.snapnum,keysel=['Group_R_Crit200','GroupFirstSub'],subcat=False)
        print "got cat"

        self.sub_ids = cat.GroupFirstSub[np.int32(self.grp_ids)]
        galf = h5py.File(self.snapbase+"/postprocessing/galprop/galprop_{}.hdf5".format(self.snapnum),'r')
        gal_sm = AU.PhysicalMass(np.array(galf['stellar_totmass']))
        x = gal_sm[self.sub_ids]

        Rvir = AU.PhysicalPosition(cat.Group_R_Crit200[np.int32(self.grp_ids)],0.19728)
        plt.scatter(np.log10(x),Rvir,marker='.',s=10,alpha=0.3,zorder=1)
        # plt.vline(150.)

        # ax.set_xscale('log')
        # ax.set_yscale('log')
        if savename != None: 
            plt.savefig(self.fig_base+savename+".pdf", bbox_inches='tight')
        else:
            plt.savefig(self.fig_base+"sm_vs_Rvir.pdf", bbox_inches='tight')

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



    def coldens_plot(self,species,kpc_mode=False,Rvir_mode=False):
        # Load all grid files.
        # Calculate radii for each one.
        # hexbin the full radius vs column density spread, by kpc or Rvir.

    
        f_npz = self.npz_base + "allgridcells_{}.npz".format(species)
        if os.path.isfile(f_npz):
            dat = np.load(f_npz)
            print list(dat)
            N = dat['N']
            r_kpc = dat['r_kpc']
            r_Rvir = dat['r_Rvir']
        else:
            i = 0
            N = np.zeros(0)
            r = np.zeros(0)
            for fn in glob.glob(self.grid_base+"s{}/*.hdf5".format(self.snapnum)):
                print "fn ",fn
                print "i ",i
                f = h5py.File(fn,'r')
                grp_id = f['Header'].attrs['grp_id']
                grp_Rvir = f['Header'].attrs['grp_Rvir']
                grid = f['grids'][species]
                grid_rad = f['Header'].attrs['grid_radius_pkpc']
                ngrid = f['Header'].attrs['ngrid']
                grid = np.array(f['grids'][species])
                f.close()

                [gridx,gridy] = np.meshgrid(np.arange(ngrid),np.arange(ngrid))
                grid_cent = (ngrid-1)/2. #assume square grid: grid_centx = grid_centy = grid_cent
                r_grid = np.sqrt((gridx-grid_cent)**2+(gridy-grid_cent)**2)
                r_kpc = self._grid_to_kpc(r_grid,ngrid,grid_rad)
                r_Rvir = r_kpc/AU.PhysicalPosition(grp_Rvir,0.19728)
                if Rvir_mode:
                    r = np.append(r,r_Rvir)
                elif kpc_mode:
                    r = np.append(r,r_kpc)
                N = np.append(N,grid)
                i+=1

            np.savez(f_npz,N=N,r_kpc=r_kpc,r_Rvir=r_Rvir)


        plt.hexbin(r_kpc,N,cmap=plt.cm.cubehelix,mincnt=1)
        plt.savefig("OVI_coldens_test.pdf")


    def plot_grids(self,species,vmin=10.,vmax=25.):

        for fn in glob.glob(self.grid_base+"s{}/*.hdf5".format(self.snapnum)):
            print "fn ",fn
            f = h5py.File(fn,'r')
            grp_id = f['Header'].attrs['grp_id']
            grid_rad = f['Header'].attrs['grid_radius_pkpc']
            ngrid = f['Header'].attrs['ngrid']
            grid = np.array(f['grids'][species])
            f.close()

            img_savepath = "{}/grids/s{}/{}_{}.pdf".format(self.fig_base,self.snapnum,str(int(grp_id)).zfill(5),species)

            maxdist = grid_rad

            plt.close('all')
            plt.imshow(grid,origin='lower',extent=(-maxdist,maxdist,-maxdist,maxdist),vmin=vmin,vmax=vmax,cmap=plt.cm.cubehelix) #spb_jet2
            bar=plt.colorbar()
            bar_label = r"log$_{10}$ N$_\mathrm{"+species+"}$ (cm$^{-2}$)"
            bar.set_label(bar_label)
            plt.xlabel(r"y (pkpc)")
            plt.ylabel(r"z (pkpc)")
            plt.savefig(img_savepath)



    def load_CGMsnap_ids(self):
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

if __name__ == '__main__':
    illustris_fan()
