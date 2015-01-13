# -*- coding: utf-8 -*-
"""Module for creating surface-density plots of various ions.  Most of this comes directly from Simeon's halohi script

Classes:
    Halo_Grid - Creates a grid around the halo center with the ion fraction calculated at each grid cell
"""

import numpy as np
import h5py
import time

import hdfsim
import cold_gas
import convert_cloudy as cc
import fieldize
#import hsml
from units import AREPO_units
AU = AREPO_units()

from mpi4py import MPI
import os.path
import readsubfHDF5


class grid_mass:

    def __init__(self):
        # self.run_list = ['']
        # self.run = 'c0_fullmetalwinds' #'gam_50_fixv'#'gam_50_fixv_nothermal' #'gam_50_BH'
        # self.base = "/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/Runs/{}/output/".format(self.run)
        # self.res=256
        # self.run = 'c0_256'#'c0_sw_256'
        # self.run_list = ['c0_256']
        # self.base_list = ["/n/hernquistfs1/Illustris/SmallBox/GFM/Production/Cosmo/Cosmo0_V6/L25n256/output/"]
        self.run_list = ['multi_speed']
        self.base_list = ["/n/home01/ptorrey/Share/apillepich/L25n256/multi_speed_winds/output/"]
        # self.run_list = ['c0_256','c2_256','c0_sw_256']
        # self.base_list = [\
        # "/n/hernquistfs1/Illustris/SmallBox/GFM/Production/Cosmo/Cosmo0_V6/L25n256/output/",\
        # "/n/hernquistfs1/Illustris/SmallBox/GFM/Production/Cosmo/Cosmo2_V6/L25n256/output/",\
        # "/n/hernquistfs1/Illustris/SmallBox/GFM/Production/Cosmo/Cosmo0_V6_strongWinds/L25n256/output/"]
        # self.run_list = ['gam_50_fixv','c0_fullmetalwinds','c0_nometalwinds','gam_50_fixv_nothermal','gam_50_BH']
        # self.base_list = [\
        # "/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/Runs/{}/output/".format(self.run_list[0]),\
        # "/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/Runs/{}/output/".format(self.run_list[1]),\
        # "/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/Runs/{}/output/".format(self.run_list[2]),\
        # "/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/Runs/{}/output/".format(self.run_list[3]),\
        # "/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/Runs/{}/output/".format(self.run_list[4]),\
        # ]
        self.res=256
        # self.snapnum_arr = np.array([60,63,68])
        # self.snapnum_arr = np.array([3,4,5])
        # self.snapnum_arr = np.array([4])
        self.snapnum_arr = np.array([15])

        self.group_min_mass = 10.**11.8 #11 CAREFUL
        self.grid_radius_pkpc = 300
        self.elem_list = ["C","C"]
        self.ion_list = [-1,4]
        # self.elem_list = ["H","Mg","C","C","Si","Si","O"]
        # self.ion_list = [1,2,3,4,3,4,6]
        # self.cloudy_type = "ion_out_fancy_atten"
        self.loadbase = '/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/CGM_new/data/CGM_snaps/'
        self.savebase = '/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/CGM_new/data/grids/'
        # self.loadbase = '/n/home04/jsuresh/'
        # self.savebase = '/n/home04/jsuresh/'

        # Here's where the magic happens
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        print "my rank = {}".format(rank)
        size = comm.Get_size()
        print "my size = {}".format(size)
        #print "done with MPI comm/rank/size initialization!"

        for rr in xrange(len(self.run_list)):
            self.run = self.run_list[rr]
            self.base = self.base_list[rr]
            for snapnum in self.snapnum_arr:
                if rank == 0:
                    # Create necessary folder if it does not exist:
                    grid_dir = self.savebase+"{}/s{}/".format(self.run,snapnum)
                    if not os.path.isdir(grid_dir):
                        print "Trying to create {}".format(grid_dir)
                        os.mkdir(grid_dir)

                    # Get header information before proceeding:
                    fname = self.base+"snapdir_"+str(snapnum).zfill(3)+"/snap_"+str(snapnum).zfill(3)+".0.hdf5"
                    f = h5py.File(fname,'r')
                    self.redshift=f["Header"].attrs["Redshift"]
                    self.hubble=f["Header"].attrs["HubbleParam"]
                    self.box=f["Header"].attrs["BoxSize"]
                    self.omegam=f["Header"].attrs["Omega0"]
                    self.omegal=f["Header"].attrs["OmegaLambda"]
                    print "redshift: ", self.redshift
                    f.close()


                    cat=readsubfHDF5.subfind_catalog(self.base,snapnum,subcat=False)
                    #cat=readsubfHDF5.subfind_catalog(base,135,keysel=["GroupMass","GroupPos","Group_R_Crit200"])

                    # Select by minimum group mass threshold:
                    self.grp_mass = np.array(cat.Group_M_Crit200) #check what hubble is
                    self.m = AU.PhysicalMass(self.grp_mass)
                    mass_select = (self.m > self.group_min_mass)

                    self.m = self.m[mass_select]
                    self.grp_mass = self.grp_mass[mass_select]
                    self.grp_ids = np.arange(np.float64(cat.ngroups))
                    self.grp_ids = self.grp_ids[mass_select]
                    self.grp_pos = np.array(cat.GroupPos)[mass_select]
                    self.grp_Rvir = np.array(cat.Group_R_Crit200)[mass_select]
                    #grp_BHmass = np.array(cat.GroupBHMass)[mass_select]
                    #grp_BHMdot = np.array(cat.GroupBHMdot)[mass_select]
                    self.n_selected_groups = np.float32(np.size(self.grp_mass))
                    del cat

                    # HERE
                    # self.grp_mass = self.grp_mass[3]
                    # self.grp_ids = self.grp_ids[3]
                    # self.grp_pos = self.grp_pos[3]
                    # self.grp_Rvir = self.grp_Rvir[3]
                    # self.n_selected_groups = 1
                else:
                    self.redshift = None
                    self.hubble = None
                    self.box = None
                    self.omegam = None
                    self.omegal = None
                    self.grp_mass = None
                    self.m = None
                    self.grp_ids = None
                    self.grp_pos = None
                    self.grp_Rvir = None
                    self.n_selected_groups = None

                # Broadcast necessary data from root process to other processes
                self.redshift = comm.bcast(self.redshift,root=0)
                self.hubble = comm.bcast(self.hubble,root=0)
                self.box = comm.bcast(self.box,root=0)
                self.omegam = comm.bcast(self.omegam,root=0)
                self.omegal = comm.bcast(self.omegal,root=0)
                self.grp_mass = comm.bcast(self.grp_mass,root=0)
                self.grp_pos = comm.bcast(self.grp_pos,root=0)
                self.m = comm.bcast(self.m,root=0)
                self.grp_ids = comm.bcast(self.grp_ids,root=0)
                self.grp_pos = comm.bcast(self.grp_pos,root=0)
                self.grp_Rvir = comm.bcast(self.grp_Rvir,root=0)
                self.n_selected_groups = comm.bcast(self.n_selected_groups,root=0)

                # Other miscellaneous stuff
                self.grid_radius = AU.CodePosition(self.grid_radius_pkpc,self.redshift,hubble=self.hubble)
                self.npart = self.res**3.
                # self.box = AU.CodePosition(self.box,self.redshift,hubble=self.hubble)  #self.box is already in code units
                self.ngrid=int(np.ceil(40*self.npart**(1./3)/self.box*2*self.grid_radius))
                # self.tab = cc.CloudyTable(self.redshift,atten_type=self.cloudy_type)
                self.tab = cc.CloudyTable(self.redshift,directory='/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/CGM/sbird/cloudy/ion_out_fancy_atten/')
                self.n_species = len(self.elem_list)


                for i in np.arange(self.n_selected_groups):
                    if (i % size == rank):
                        print "Projecting i / task:", i, rank
                        grp_id = self.grp_ids[i]
                        gpos  = self.grp_pos[i]
                        r200  = self.grp_Rvir[i]
                        
                        self.grid = np.zeros([self.n_species,self.ngrid,self.ngrid])
                        #grid = np.zeros([ngrid,ngrid])

                        CGMsnap_file_path = self.loadbase+ "{}/s{}/{}.hdf5".format(self.run,snapnum,str(int(grp_id)).zfill(5))
                        outputpath = self.savebase+ "{}/s{}/".format(self.run,snapnum)
                        savename = outputpath+str(int(grp_id)).zfill(5)+'.hdf5'

                        if not os.path.isfile(CGMsnap_file_path):
                            print "File {} does not exist!  Skipping...".format(CGMsnap_file_path)
                            print "[Note: log10(mass) = {}]".format(np.log10(self.m[i]))
                        # elif os.path.isfile(savename):
                        #     print "File {} already exists!  Skipping...".format(savename)
                        else:
                            print "working on {}".format(savename)
                            [m,pos,metals,rho,u,nelec,hsml,T,neut_frac] = self.load_CGM_file(CGMsnap_file_path,i)

                            pos_cent = self._fix_pos(self.grp_pos[i],pos,self.box)
                            self.calc_grid(m,pos_cent,metals,rho,hsml,T,neut_frac,kernel_type='SPH')

                            # print "Not saving! "
                            if True:
                                print "Saving group file now to {}".format(savename)
                                print ""
                                # do not overwrite mode:
                                f=h5py.File(savename,'w-')

                                grp = f.create_group("Header")
                                grp.attrs["hubble"]=self.hubble
                                grp.attrs["omegam"]=self.omegam
                                grp.attrs["omegal"]=self.omegal
                                grp.attrs["redshift"]=self.redshift
                                grp.attrs["box"]=self.box

                                grp.attrs["grp_id"] = grp_id
                                grp.attrs["grp_mass"] = self.grp_mass[i]
                                grp.attrs["grp_pos"] = self.grp_pos[i]
                                grp.attrs["grp_Rvir"] = self.grp_Rvir[i]
                                #grp.attrs["grp_BHmass"] = grp_BHmass[i]
                                #grp.attrs["grp_BHMdot"] = grp_BHMdot[i]

                                grp.attrs["grid_radius_pkpc"] = self.grid_radius_pkpc
                                grp.attrs["ngrid"] = self.ngrid

                                p_grp = f.create_group('grids')
                                for i in xrange(self.n_species):
                                    elem = self.elem_list[i]
                                    ion = self.ion_list[i]
                                    dataset_name = elem+str(ion)
                                    p_grp.create_dataset(dataset_name,data=self.grid[i,:,:])
                                f.close()




    def load_CGM_file(self,CGMsnap_file_path,i,use_block_name=False):
        f=h5py.File(CGMsnap_file_path,'r')
        print "CGMsnap_file_path: ",CGMsnap_file_path
        bar = f['PartType0']

        #print list(bar)
        #dat_str_list = ["Masses","Coordinates","GFM_Metals","GFM_Metallicity","Velocities","Density","Volume","InternalEnergy","ElectronAbundance","Radius"]
        m = np.array(bar["Masses"])
        pos = np.array(bar["Coordinates"])
        metals = np.array(bar["GFM_Metals"])
        rho = np.array(bar["Density"])
        # hsml = np.array(bar["SmoothingLength"])
        u = np.array(bar["InternalEnergy"])
        nelec = np.array(bar["ElectronAbundance"])
        neut_frac = np.array(bar["NeutralHydrogenAbundance"])

        vol = np.array(bar["Volume"])
        hsml = (3.*vol/4./np.pi)**(0.33333333)
        hsml *= 2. # SPH fudge factor [see research notes from 10/10/14]
        del vol
        T = AU.GetTemp(u, nelec, gamma = 5.0/3.0)

        # Modify positions if they are on the other side of the box:
        for j in np.arange(3):
            bc1 = np.abs(pos[:,j]-self.grp_pos[i,j]-self.box) < np.abs(pos[:,j]-self.grp_pos[i,j])
            pos[bc1,j] -= self.box

            bc2 = np.abs(pos[:,j]-self.grp_pos[i,j]+self.box) < np.abs(pos[:,j]-self.grp_pos[i,j])
            pos[bc2,j] += self.box

        return [m,pos,metals,rho,u,nelec,hsml,T,neut_frac]


    def sub_gridize_single_halo(self,pos,hsml,ion_mass,grid,weights=None,kernel_type='SPH'):
        coords=fieldize.convert_centered(pos,self.ngrid,2*self.grid_radius)

        #Convert smoothing lengths to grid coordinates.
        hsml_grid = hsml*(self.ngrid/(2*self.grid_radius))
        
        #interpolate the density
        ts = time.time()
        print "starting to fieldize"
        if kernel_type == 'SPH':
            fieldize.sph_str(coords,ion_mass,grid,hsml_grid,weights=weights)
        elif kernel_type == 'TOPHAT':
            fieldize.tophat_str(coords,ion_mass,grid,hsml_grid,weights=weights)
        print "time to fieldize: ",time.time()-ts
        return


    def calc_grid(self,m,pos,metals,rho,hsml,T,neut_frac,kernel_type='SPH'):
        print "initial # of particles before picking out: ",np.size(m)
        for d in np.arange(3):
            ind = np.abs(pos[:,d]) <= self.grid_radius
            pos = pos[ind]
            m = m[ind]
            metals = metals[ind]
            hsml = hsml[ind]
            T = T[ind]
            rho = rho[ind]  
            neut_frac = neut_frac[ind]
            #print "np.sum(ind) ",np.sum(ind)

        print "# of particles to fieldize over: ",np.size(m)

        # Impose cut on T and rho so cloudy doesn't die:
        T_cut = np.logical_and(np.log10(T) >= 3.,np.log10(T) <= 8.6)
        rho_phys = AU.PhysicalDensity(rho,self.redshift,hubble=self.hubble)
        rho_Hatoms = rho_phys*(metals[:,0]/AU.ProtonMass)
        rho_cut = np.logical_and(np.log10(rho_Hatoms) >= -7.,np.log10(rho_Hatoms) <= 4.)
        cloudy_cut = np.logical_and(T_cut,rho_cut)

        if np.sum(cloudy_cut) < np.size(m): 
            print "cloudy cut removed {} particles".format(np.size(m)-np.sum(cloudy_cut))
            pos = pos[cloudy_cut]
            m = m[cloudy_cut]
            metals = metals[cloudy_cut]
            hsml = hsml[cloudy_cut]
            T = T[cloudy_cut]
            rho = rho[cloudy_cut]
            rho_Hatoms = rho_Hatoms[cloudy_cut]
            neut_frac = neut_frac[cloudy_cut] 

        for i in xrange(self.n_species):
            elem = self.elem_list[i]
            ion = self.ion_list[i]
            print "species: {}".format(elem+str(ion))


            # elem_massfrac = metals[:,AU.elem_lookup(elem)]
            # ion_frac = self.tab.ion(elem,ion,rho_Hatoms,T)
            # species_mass = m*elem_massfrac*ion_frac

            if ion == -1:
                # Take full element, e.g. all carbon
                elem_massfrac = metals[:,AU.elem_lookup(elem)]
                species_mass = m*elem_massfrac

            elif elem == "H" and ion == 1:
                H_massfrac = metals[:,0]
                star=cold_gas.RahmatiRT(self.redshift, self.hubble)
                fake_bar = {'Density':rho,'NeutralHydrogenAbundance':neut_frac}
                new_neut_frac = star.get_reproc_HI(fake_bar)
                species_mass = m*H_massfrac * new_neut_frac

            else: 
                elem_massfrac = metals[:,AU.elem_lookup(elem)]
                ion_frac = self.tab.ion(elem,ion,rho_Hatoms,T)
                species_mass = m*elem_massfrac*ion_frac

            self.sub_gridize_single_halo(pos,hsml,species_mass,self.grid[i,:,:],kernel_type=kernel_type)
            self.grid[i,:,:][self.grid[i,:,:] < 0.] = 0.


            elem_mass_g = AU.elem_atom_mass(elem)
            massg=AU.UnitMass_in_g/self.hubble*(1/elem_mass_g)
            epsilon=2.*self.grid_radius/self.ngrid*AU.UnitLength_in_cm/self.hubble/(1+self.redshift)
            self.grid[i,:,:] *= (massg/epsilon**2)
            self.grid[i,:,:] += 0.1
            np.log10(self.grid[i,:,:],self.grid[i,:,:])


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



            

if __name__ == '__main__':
    grid_mass()
    
# if __name__ == '__main__':
#     main()








