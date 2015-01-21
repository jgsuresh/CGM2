# -*- coding: utf-8 -*-
"""Module for creating surface-density plots of various ions.  
Most of this comes directly from Simeon's halohi script

Classes:
    species_grid - Creates a grid around the halo center with the ion fraction calculated at each grid cell
"""

import numpy as np
import h5py
import time
from mpi4py import MPI
import os.path

import readsubfHDF5
import hdfsim
import cold_gas
import convert_cloudy as ccl
import fieldize

from units import AREPO_units
AU = AREPO_units()

class species_grid:
    def __init__(self,run=None,snap_base=None,CGMsnap_base=None,save_base=None,res=None,snapnum=None,grp_ids=None,\
        grid_radius_pkpc=None,elem_list=None,ion_list=None,cloudy_type=None,cloudy_dir=None,\
        verbose=False,**kwargs):
        self.run = run
        self.snap_base = snap_base
        self.CGMsnap_base = CGMsnap_base
        self.save_base = save_base
        self.res = res
        self.snapnum = snapnum
        # self.group_min_mass = group_min_mass
        self.grid_radius_pkpc = grid_radius_pkpc
        self.elem_list = elem_list
        self.ion_list = ion_list
        self.cloudy_type = cloudy_type
        self.cloudy_dir = cloudy_dir
        self.kwargs = kwargs
        self.grp_ids = grp_ids
        self.verbose = verbose

        # Or hardcode values:
        # self.run = None
        # self.snap_base = None
        # self.CGMsnap_base = '/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/CGM_new/data/CGM_snaps/'
        # self.save_base = '/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/CGM_new/data/grids/'
        # self.res = None #Future: get this from the simulation
        # self.snapnum = None
        # # self.group_min_mass = None
        # self.grid_radius_pkpc = None
        # self.elem_list = ["C","C"]
        # self.ion_list = [-1,4]
        # self.cloudy = "UVB_sf_xrays_ext"
        # self.cloudy_dir = '/n/home04/jsuresh/scratch1/Cloudy_test/UVB_sf_xrays_ext/'

        # Here's where the magic happens
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        if self.verbose:
            print "my rank = {}".format(rank)
            print "my size = {}".format(size)
        #print "done with MPI comm/rank/size initialization!"

        # Create necessary folder if it does not exist:
        if rank == 0:
            grid_dir = self.save_base+"{}/s{}/".format(self.run,snapnum)
            if not os.path.isdir(grid_dir):
                if self.verbose:
                    print "Creating {}".format(grid_dir)
                os.mkdir(grid_dir)

            # Load header information from snapshot:
            self.load_header()
        else:
            self.redshift = None
            self.hubble = None
            self.box = None

        # Broadcast necessary data from root process to other processes
        self.redshift = comm.bcast(self.redshift,root=0)
        self.hubble = comm.bcast(self.hubble,root=0)
        self.box = comm.bcast(self.box,root=0)

        # Other set-up
        self.grid_radius = AU.CodePosition(self.grid_radius_pkpc,self.redshift,hubble=self.hubble)
        self.npart = self.res**3.
        self.ngrid=int(np.ceil(40*self.npart**(1./3)/self.box*2*self.grid_radius))
        tab = ccl.CloudyTable(self.redshift,directory=cloudy_dir)
        self.n_species = len(self.elem_list)
        self.n_selected_groups = np.size(self.grp_ids)

        if self.cloudy_type == "UVB_sf_xrays_ext":
            if not kwargs.has_key(multiple_sources):
                try:
                    self.gal_SFR = kwargs['gal_SFR']
                except:
                    raise KeyError("Need gal_SFR to be passed in for cloudy type {}".format(self.cloudy_type))

            if kwargs.has_key(multiple_sources) and kwargs['multiple_sources'] == True:
                try:
                    self.sub_id_dict = kwargs['sub_id_dict'] # Dictionary which, for each group, has list of corresponding subhalo IDs
                    self.sub_pos = kwargs['sub_pos']
                    self.sub_SM = kwargs['sub_SM']
                    self.sub_SFR = kwargs['sub_SFR']
                    self.grp_firstsub = kwargs['grp_firstsub']

                except:
                    raise KeyError("Missing some of the subhalo data which are necessary to implement multiple sources.")


        for i in np.arange(self.n_selected_groups):
            if (i % size == rank):
                if self.verbose:
                    print "Projecting i / task:", i, rank
                grp_id = self.grp_ids[i]

                # Work on this group if CGM snapshot exists and grid file does not already exist:
                CGMsnap_file_path = self.CGMsnap_base + "{}/s{}/{}.hdf5".format(self.run,self.snapnum,str(int(grp_id)).zfill(5))
                save_path = self.save_base + "{}/s{}/{}.hdf5".format(self.run,snapnum,str(int(grp_id)).zfill(5))
                if not os.path.isfile(CGMsnap_file_path):
                    if self.verbose: 
                        print "File {} does not exist!  Skipping...".format(CGMsnap_file_path)
                elif os.path.isfile(save_path):
                    if self.verbose: 
                        print "File {} already exists!  Skipping...".format(save_path)
                else:
                    if self.verbose: 
                        print "Working on {}".format(save_path)
                    data_dict = self.load_CGM_file(grp_id,CGMsnap_file_path,center_positions=True,radial_cut=True)
                
                    grid_dict = self.calc_grid(data_dict,kernel_type='SPH')
                    self.save_grid(grp_id,grid_dict,save_path)

                    

    #==================================================================================================
    def load_header(self):
        fname = self.snap_base+"snapdir_"+str(self.snapnum).zfill(3)+"/snap_"+str(self.snapnum).zfill(3)+".0.hdf5"
        f = h5py.File(fname,'r')
        self.redshift=f["Header"].attrs["Redshift"]
        self.hubble=f["Header"].attrs["HubbleParam"]
        self.box=f["Header"].attrs["BoxSize"]
        self.omegam=f["Header"].attrs["Omega0"]
        self.omegal=f["Header"].attrs["OmegaLambda"]
        print "redshift: ", self.redshift
        f.close()

    #==================================================================================================
    def load_CGM_file(self,grp_id,CGMsnap_file_path,use_block_name=False,center_positions=True,radial_cut=True,cloudy_cut=True):
        f=h5py.File(CGMsnap_file_path,'r')
        
        grp_pos = f['Header'].attrs["grp_pos"]
        grp_Rvir = f['Header'].attrs["grp_Rvir"]

        bar = f['PartType0']
        # dat_str_list = ["Masses","Coordinates","GFM_Metals","GFM_Metallicity","Velocities","Density","Volume","InternalEnergy","ElectronAbundance","NeutralHydrogenAbundance","SmoothingLength"]
        # self.block_list = ["MASS","POS ","GMET","GZ  ","VEL ","RHO ","VOL ","U   ","NE  ","NH  ","HSML"]
        if use_block_name:
            m = np.array(bar["MASS"])
            pos = np.array(bar["POS "])
            metals = np.array(bar["GMET"])
            rho = np.array(bar["RHO "])
            u = np.array(bar["U   "])
            nelec = np.array(bar["NE  "])
            neut_frac = np.array(bar["NH  "])
            vol = np.array(bar["VOL "])
        else:
            m = np.array(bar["Masses"])
            pos = np.array(bar["Coordinates"])
            metals = np.array(bar["GFM_Metals"])
            rho = np.array(bar["Density"])
            # hsml = np.array(bar["SmoothingLength"])
            u = np.array(bar["InternalEnergy"])
            nelec = np.array(bar["ElectronAbundance"])
            neut_frac = np.array(bar["NeutralHydrogenAbundance"])
            vol = np.array(bar["Volume"])
        f.close()

        hsml = (3.*vol/4./np.pi)**(0.33333333)
        hsml *= 2. # SPH fudge factor [see research notes from 10/10/14]
        del vol
        T = AU.GetTemp(u, nelec, gamma = 5.0/3.0)

        # Modify positions if they are on the other side of the box:
        for j in np.arange(3):
            bc1 = np.abs(pos[:,j]-grp_pos[j]-self.box) < np.abs(pos[:,j]-grp_pos[j])
            pos[bc1,j] -= self.box

            bc2 = np.abs(pos[:,j]-grp_pos[j]+self.box) < np.abs(pos[:,j]-grp_pos[j])
            pos[bc2,j] += self.box

        # Center positions on group position
        if center_positions:
            pos_cent = self._fix_pos(grp_pos,pos,self.box)
            pos = pos_cent

        # Implement radial cut for particles on grid
        if radial_cut:
            if self.verbose: print "Initial # of particles before doing radial cut: ",np.size(m)

            for d in np.arange(3):
                ind = np.abs(pos[:,d]) <= 1.5*self.grid_radius
                pos = pos[ind]
                m = m[ind]
                metals = metals[ind]
                hsml = hsml[ind]
                T = T[ind]
                rho = rho[ind]  
                neut_frac = neut_frac[ind]

            if self.verbose: print "Final # of particles after doing radial cut: ",np.size(m)

        # Impose cut on T and rho so cloudy doesn't die:
        rho_phys = AU.PhysicalDensity(rho,self.redshift,hubble=self.hubble)
        rho_Hatoms = rho_phys*(metals[:,0]/AU.ProtonMass)
        if cloudy_cut:
            T_cut = np.logical_and(np.log10(T) >= 3.,np.log10(T) <= 8.6)
            rho_cut = np.logical_and(np.log10(rho_Hatoms) >= -7.,np.log10(rho_Hatoms) <= 4.)
            cloudy_cut = np.logical_and(T_cut,rho_cut)

            if np.sum(cloudy_cut) < np.size(m): 
                if self.verbose: print "Cloudy cut removed {} particles".format(np.size(m)-np.sum(cloudy_cut))
                pos = pos[cloudy_cut]
                m = m[cloudy_cut]
                metals = metals[cloudy_cut]
                hsml = hsml[cloudy_cut]
                T = T[cloudy_cut]
                rho = rho[cloudy_cut]
                rho_Hatoms = rho_Hatoms[cloudy_cut]
                neut_frac = neut_frac[cloudy_cut] 

        # return [m,pos,metals,rho,u,nelec,hsml,T,neut_frac]
        data_dict = {}
        data_dict['m'] = m
        data_dict['pos'] = pos
        data_dict['metals'] = metals
        data_dict['rho'] = rho
        data_dict['rho_Hatoms'] = rho_Hatoms
        data_dict['u'] = u
        data_dict['nelec'] = nelec
        data_dict['hsml'] = hsml
        data_dict['T'] = T
        data_dict['neut_frac'] = neut_frac
        data_dict['grp_pos'] = grp_pos
        data_dict['grp_id'] = grp_id
        return data_dict


    #==================================================================================================
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


    #==================================================================================================
    def calc_grid(self,data_dict,kernel_type='SPH'):
        m = data_dict['m']
        pos = data_dict['pos']
        metals = data_dict['metals']
        rho = data_dict['rho']
        rho_Hatoms = data_dict['rho_Hatoms']
        u = data_dict['u']
        nelec = data_dict['nelec']
        hsml = data_dict['hsml']
        T = data_dict['T']
        neut_frac = data_dict['neut_frac']
        grp_id = data_dict['grp_id']
        grp_pos = data_dict['grp_pos']

        r = AU.PhysicalPosition(np.sqrt(pos[:,0]**2.+pos[:,1]**2.+pos[:,2]**2.),self.redshift)

        if self.verbose: print "# of particles to fieldize over: ",np.size(m)

        grid_dict = {}
        for ss in xrange(self.n_species):
            elem = self.elem_list[ss]
            ion = self.ion_list[ss]
            if self.verbose:
                print "Species: {}".format(elem+str(ion))
                print "Cloudy: {}".format(self.cloudy_type)

            # For neutral hydrogen, use Rahmati fitting formula to compute HI fraction
            if elem == "H" and ion == 1:
                H_massfrac = metals[:,0]
                star=cold_gas.RahmatiRT(self.redshift, self.hubble)
                fake_bar = {'Density':rho,'NeutralHydrogenAbundance':neut_frac}
                new_neut_frac = star.get_reproc_HI(fake_bar)
                species_mass = m*H_massfrac * new_neut_frac

            # Otherwise, use cloudy to compute ion fraction
            else:
                elem_massfrac = metals[:,AU.elem_lookup(elem)]

                if self.cloudy_type == "ion_out_fancy_atten":
                    ion_frac = self.tab.ion(elem,ion,rho_Hatoms,T)
                elif self.cloudy_type == "UVB_sf_xrays_ext":

                    if self.kwargs.has_key('multiple_sources') and self.kwargs['multiple_sources'] == True:
                        # Count # of subhalos for this group:
                        sub_ids = self.sub_id_dict[grp_id]
                        grp_firstsubid = self.grp_firstsub[grp_id]
                        n_sources = np.size(sub_ids)

                        if n_sources == 0:
                            raise ValueError("No subhalos for group {}!  ERROR".format(grp_id))
                        elif n_sources == 1:
                            beta = self.sub_SFR[grp_firstsubid]/r**2.
                        else:
                            # loop over all sources; for each source, get its position, SFR, and stellar mass
                            beta = np.zeros_like(r)
                            # beta_arr = np.zeros([n_sources,np.size(r)])
                            for sub_id in sub_ids:
                                sub_SFR = self.sub_SFR[sub_id]
                                sub_SM = self.sub_SM[sub_id]
                                sub_pos = self.sub_pos[sub_id]

                                if sub_SFR == 0: pass
                                else:
                                    # Put subhalo position in frame where group center is at origin
                                    subpos_cent = self._fix_pos(grp_pos,np.array([sub_pos]),self.box)
                                    subpos_cent = subpos_cent[0]
                                    # Then calculate distance from each cell to the subhalo position
                                    rs = AU.PhysicalPosition(np.sqrt(\
                                    (pos[:,0]-subpos_cent[0])**2.+\
                                    (pos[:,1]-subpos_cent[1])**2.+\
                                    (pos[:,2]-subpos_cent[2])**2.)\
                                    ,self.redshift)

                                    # Sanity check:
                                    if np.sum(rs > self.box/2.) > 0: raise ValueError("Radius error!")

                                    beta += sub_SFR/rs**2.
                            
                    else:
                        beta = self.gal_SFR[grp_id]/r**2.

                    if self.verbose:
                        print "beta {}".format(beta)
                    # Given beta for each cell (radiation that each cell sees from young stellar pops), calculate 
                    # ionization states
                    ion_frac = tab.ion(elem,ion,rho_Hatoms,T,beta=beta)

                # Given element mass fraction and ion fraction within that element, we can calculate the species mass in each cell:
                species_mass = m*elem_massfrac*ion_frac

            grid = np.zeros([self.ngrid,self.ngrid])
            self.sub_gridize_single_halo(pos,hsml,species_mass,grid,kernel_type=kernel_type)

            # Put grid in correct units
            elem_mass_g = AU.elem_atom_mass(elem)
            massg=AU.UnitMass_in_g/self.hubble*(1/elem_mass_g)
            epsilon=2.*self.grid_radius/self.ngrid*AU.UnitLength_in_cm/self.hubble/(1+self.redshift)
            grid *= (massg/epsilon**2)
            grid += 0.1
            np.log10(grid,grid)

            grid_dict[elem+str(ion)] = np.copy(grid)

        return grid_dict


    #==================================================================================================
    def sub_gridize_single_halo(self,pos,hsml,ion_mass,grid,weights=None,kernel_type='SPH'):
        coords=fieldize.convert_centered(pos,self.ngrid,2*self.grid_radius)

        #Convert smoothing lengths to grid coordinates.
        hsml_grid = hsml*(self.ngrid/(2*self.grid_radius))
        
        #interpolate the density
        ts = time.time()
        if self.verbose:
            print "Starting to fieldize"
        if kernel_type == 'SPH':
            fieldize.sph_str(coords,ion_mass,grid,hsml_grid,weights=weights)
        elif kernel_type == 'TOPHAT':
            fieldize.tophat_str(coords,ion_mass,grid,hsml_grid,weights=weights)
        if self.verbose:
            print "Time to fieldize: ",time.time()-ts
        return


    #==================================================================================================
    def save_grid(self,grp_id,grid_dict,save_path):
        if self.verbose: 
            print "Saving group file now to {}".format(save_path)

        f=h5py.File(save_path,'w-')

        grp = f.create_group("Header")
        grp.attrs["hubble"]=self.hubble
        grp.attrs["omegam"]=self.omegam
        grp.attrs["omegal"]=self.omegal
        grp.attrs["redshift"]=self.redshift
        grp.attrs["box"]=self.box

        grp.attrs["grp_id"] = grp_id
        grp.attrs["grid_radius_pkpc"] = self.grid_radius_pkpc
        grp.attrs["ngrid"] = self.ngrid
        grp.attrs["cloudy_type"] = self.cloudy_type

        grids = f.create_group("grids")
        for key in grid_dict:
            grids.create_dataset(key,data=grid_dict[key])
        f.close()
