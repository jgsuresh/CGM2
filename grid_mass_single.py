# -*- coding: utf-8 -*-
"""Module for creating surface-density plots of various ions.  Most of this comes directly from Simeon's halohi script

Classes:
    Halo_Grid - Creates a grid around the halo center with the ion fraction calculated at each grid cell
"""

import numpy as np
import h5py
import time

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
        self.res = 1820
        self.box = 75000
        self.redshift = 0.2
        self.hubble = 0.7
        fn = "/n/home04/jsuresh/scratch1/Cloudy_test/cutout/CGM_snap/02774_full300.hdf5"
        # fn = "/n/home04/jsuresh/scratch1/Cloudy_test/cutout/02774.hdf5"
        # fn = "/n/ghernquist/ptorrey/Illustris/GroupParsedSnapshots/Illustris-1/snapshot_120/subfolder_027/group_2774.hdf5"
        
        self.grid_radius_pkpc = 300
        # self.elem_list = ["H","C","C","C","N","N","N","N","O","O","O","O","O","O","Ne","Mg","Mg","Si","Si","Si","Fe"]
        # self.ion_list = [1,2,3,4,2,3,4,5,3,4,5,6,7,8,8,2,10,3,4,12,2]
        self.elem_list = ["H","C","O"]
        self.ion_list = [1,4,6]
        # cloudy_dir
        # self.cloudy_type = "ion_out_fancy_atten"
        self.cloudy_type = "UVB_sf_xrays_ext"
        self.SFR = 100. # Msun/yr
        # self.beta = 10.**-4.
        self.savebase = '/n/home04/jsuresh/scratch1/Cloudy_test/cutout/grids/full_cutout/'

        # Other miscellaneous stuff
        self.grid_radius = AU.CodePosition(self.grid_radius_pkpc,self.redshift,hubble=self.hubble)
        self.npart = self.res**3.
        self.box = AU.CodePosition(self.box,self.redshift,hubble=self.hubble) #HERE
        self.ngrid=int(np.ceil(40*self.npart**(1./3)/self.box*2*self.grid_radius))
        # self.tab = cc.CloudyTable(self.redshift,atten_type=self.cloudy_type)
        self.tab = cc.CloudyTable(0.2,directory='/n/home04/jsuresh/scratch1/Cloudy_test/UVB_sf_xrays_ext/')
        self.n_species = len(self.elem_list)

        self.grid = np.zeros([self.n_species,self.ngrid,self.ngrid])
        f = h5py.File(fn,'r')
        self.grp_pos = f['Header'].attrs['grp_pos']
        f.close()
        [m,pos,metals,rho,u,nelec,hsml,T,neut_frac] = self.load_CGM_file(fn,use_block_name=False)
        

        pos_cent = self._fix_pos(self.grp_pos,pos,self.box)
        print "pos_cent: ",pos_cent
        self.calc_grid(m,pos_cent,metals,rho,hsml,T,neut_frac,kernel_type='SPH')

        # savename = self.savebase + self.cloudy_type + "_beta{}.hdf5".format(self.beta)
        savename = self.savebase + "g{}_".format(2774)+self.cloudy_type + "_interp_SFR{}.hdf5".format(self.SFR)
        print "Saving group file now to {}".format(savename)
        f=h5py.File(savename,'w-')
        grp = f.create_group("Header")
        grp.attrs["grid_radius_pkpc"] = self.grid_radius_pkpc
        grp.attrs["ngrid"] = self.ngrid

        p_grp = f.create_group('grids')
        for i in xrange(self.n_species):
            elem = self.elem_list[i]
            ion = self.ion_list[i]
            dataset_name = elem+str(ion)
            p_grp.create_dataset(dataset_name,data=self.grid[i,:,:])
        f.close()
        




    def load_CGM_file(self,CGMsnap_file_path,use_block_name=True):
        f=h5py.File(CGMsnap_file_path,'r')
        bar = f['PartType0']

        #print list(bar)
        #dat_str_list = ["Masses","Coordinates","GFM_Metals","GFM_Metallicity","Velocities","Density","Volume","InternalEnergy","ElectronAbundance","Radius"]
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

        hsml = (3.*vol/4./np.pi)**(0.33333333)
        hsml *= 2.
        del vol
        T = AU.GetTemp(u, nelec, gamma = 5.0/3.0)

        # Modify positions if they are on the other side of the box:
        for j in np.arange(3):
            bc1 = np.abs(pos[:,j]-self.grp_pos[j]-self.box) < np.abs(pos[:,j]-self.grp_pos[j])
            pos[bc1,j] -= self.box

            bc2 = np.abs(pos[:,j]-self.grp_pos[j]+self.box) < np.abs(pos[:,j]-self.grp_pos[j])
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
        rho_phys = AU.PhysicalDensity(rho,self.redshift,hubble=self.hubble)
        rho_Hatoms = rho_phys*(metals[:,0]/AU.ProtonMass)
        rho_cut = np.logical_and(np.log10(rho_Hatoms) >= -7.,np.log10(rho_Hatoms) <= 4.)
        if False:
            low_dens = rho_Hatoms < 10.**-3
            T[low_dens] *= 10.**0.5
        T_cut = np.logical_and(np.log10(T) >= 3.,np.log10(T) <= 8.6)
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



        # Calculate radial distance for particles from center:
        r = AU.PhysicalPosition(np.sqrt(pos[:,0]**2.+pos[:,1]**2.+pos[:,2]**2.),self.redshift,hubble=self.hubble)
        beta = self.SFR/(r**2.)

        for i in xrange(self.n_species):
            elem = self.elem_list[i]
            ion = self.ion_list[i]
            print "species: {}".format(elem+str(ion))


            # elem_massfrac = metals[:,AU.elem_lookup(elem)]
            # ion_frac = self.tab.ion(elem,ion,rho_Hatoms,T)
            # species_mass = m*elem_massfrac*ion_frac

            if elem == "H" and ion == 1:
                H_massfrac = metals[:,0]
                star=cold_gas.RahmatiRT(self.redshift, self.hubble)
                fake_bar = {'Density':rho,'NeutralHydrogenAbundance':neut_frac}
                new_neut_frac = star.get_reproc_HI(fake_bar)
                species_mass = m*H_massfrac * new_neut_frac

            else: 
                elem_massfrac = metals[:,AU.elem_lookup(elem)]
                # ion_frac = self.tab.ion(elem,ion,rho_Hatoms,T,beta=self.beta*np.ones_like(T))
                ion_frac = self.tab.ion(elem,ion,rho_Hatoms,T,beta=beta)
                species_mass = m*elem_massfrac*ion_frac

            self.sub_gridize_single_halo(pos,hsml,species_mass,self.grid[i,:,:],kernel_type=kernel_type)

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








