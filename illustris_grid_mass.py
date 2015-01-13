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
import convert_cloudy as ccl
import fieldize
#import hsml
from units import AREPO_units
AU = AREPO_units()

from mpi4py import MPI
import os.path
import readsubfHDF5

import glob


class grid_mass:

    def __init__(self):
        self.base = "/n/ghernquist/Illustris/Runs/Illustris-1/output/"
        self.res=1820
        self.snapnum_arr = np.array([120])
        self.redshift = 0.19728
        self.grid_radius_pkpc = 200
        # self.elem_list = ["H","Mg","Si","N","O"]
        self.elem_list = ["N","O"]
        self.ion_list = [5,6]
        # self.cloudy_list = ["ion_out_fancy_atten","UVB_sf_xrays_ext"]
        self.cloudy_list = ["UVB_sf_xrays_ext"]
        self.tab_list = [\
        ccl.CloudyTable(self.redshift,directory='/n/home04/jsuresh/scratch1/Cloudy_test/UVB_sf_xrays_ext/')]
        # ccl.CloudyTable(self.redshift,directory='/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/CGM/sbird/cloudy/ion_out_fancy_atten/'),
        self.n_cloudy_var = len(self.cloudy_list)

        # Modify temperatures of low-density CGM?
        self.modify_temp = False
        self.tempfac_arr = np.array([1,3.,10.])
        self.n_tempfac = len(self.tempfac_arr)
        
        # Multiple sources of radiation, or only use central galaxy?
        self.multiple_sources = True
        # self.mmin_source = 10.**9.

        if self.modify_temp and self.multiple_sources:
            raise NotImplementedError

        self.loadbase = '/n/home04/jsuresh/scratch1/AREPOfest/data/CGM_snaps/'
        self.savebase = '/n/home04/jsuresh/scratch1/AREPOfest/data/grids/'
        # self.loadbase = '/n/home04/jsuresh/scratch1/QCGM2/data/CGM_snaps/'
        # self.savebase = '/n/home04/jsuresh/scratch1/QCGM2/data/grids/'
        


        # Here's where the magic happens
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        print "my rank = {}".format(rank)
        size = comm.Get_size()
        print "my size = {}".format(size)
        #print "done with MPI comm/rank/size initialization!"

        for snapnum in self.snapnum_arr:
            if rank == 0:
                # Create necessary folder if it does not exist:
                grid_dir = self.savebase+"s{}/".format(snapnum)
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

                self.grp_ids = np.zeros(0)
                self.grp_mass = np.zeros(0)
                self.grp_pos = np.zeros([0,3])
                self.grp_Rvir = np.zeros(0)
                if True:
                    for fn in glob.glob(self.loadbase+"s{}/*.hdf5".format(snapnum)):
                        # print "fn ",fn
                        f = h5py.File(fn,'r')
                        grp_id = f['Header'].attrs['grp_id']
                        m = AU.PhysicalMass(f['Header'].attrs['grp_mass'])
                        x = f['Header'].attrs['grp_pos']
                        R = f['Header'].attrs['grp_Rvir']
                        f.close()

                        x = np.array([x])
                        self.grp_ids = np.append(self.grp_ids,grp_id)
                        self.grp_mass = np.append(self.grp_mass,m)
                        self.grp_pos = np.append(self.grp_pos,x,axis=0)
                        self.grp_Rvir = np.append(self.grp_Rvir,R)
                else:
                    fn = "/n/home04/jsuresh/scratch1/AREPOfest/data/CGM_snaps/s120/00300.hdf5"
                    print "fn ",fn
                    f = h5py.File(fn,'r')
                    grp_id = f['Header'].attrs['grp_id']
                    m = AU.PhysicalMass(f['Header'].attrs['grp_mass'])
                    x = f['Header'].attrs['grp_pos']
                    R = f['Header'].attrs['grp_Rvir']
                    f.close()

                    x = np.array([x])
                    self.grp_ids = np.append(self.grp_ids,grp_id)
                    self.grp_mass = np.append(self.grp_mass,m)
                    self.grp_pos = np.append(self.grp_pos,x,axis=0)
                    self.grp_Rvir = np.append(self.grp_Rvir,R)
                

                foo = np.argsort(self.grp_ids)
                self.grp_ids = self.grp_ids[foo]
                self.grp_mass = self.grp_mass[foo]
                self.grp_pos = self.grp_pos[foo]
                self.grp_Rvir = self.grp_Rvir[foo]
                self.n_selected_groups = np.float32(np.size(self.grp_mass))

            else:
                self.redshift = None
                self.hubble = None
                self.box = None
                self.omegam = None
                self.omegal = None
                self.grp_mass = None
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
            self.grp_ids = comm.bcast(self.grp_ids,root=0)
            self.grp_pos = comm.bcast(self.grp_pos,root=0)
            self.grp_Rvir = comm.bcast(self.grp_Rvir,root=0)
            self.n_selected_groups = comm.bcast(self.n_selected_groups,root=0)

            # Other miscellaneous stuff
            self.grid_radius = AU.CodePosition(self.grid_radius_pkpc,self.redshift,hubble=self.hubble)
            self.npart = self.res**3.
            self.ngrid=int(np.ceil(40*self.npart**(1./3)/self.box*2*self.grid_radius))
            # self.tab = cc.CloudyTable(self.redshift,atten_type=self.cloudy_type)
            self.n_species = len(self.elem_list)


            if "UVB_sf_xrays_ext" in self.cloudy_list:
                galprop = h5py.File(self.base+"postprocessing/galprop/galprop_{}.hdf5".format(snapnum),'r')
                cat = readsubfHDF5.subfind_catalog(self.base,snapnum,keysel=['GroupFirstSub'])
                gal_SFR = np.array(galprop['gas_sfr_inrad'])

                cent_ids = cat.GroupFirstSub[np.int32(self.grp_ids)]
                self.grp_SFR = np.array(galprop['gas_sfr_inrad'])[cent_ids]

                if self.multiple_sources:
                    cat_sat = readsubfHDF5.subfind_catalog(self.base,snapnum,keysel=['SubhaloPos','SubhaloGrNr'])
                    nsubs = np.size(cat_sat.SubhaloGrNr)
                    sub_ids = np.arange(nsubs)

                    # get ids for all subhalos associated with each group
                    self.subhalo_id_dict = {}
                    for grnr in self.grp_ids:
                        gr_sub_ids = np.copy(sub_ids[cat_sat.SubhaloGrNr == grnr])
                        self.subhalo_id_dict[grnr] = gr_sub_ids

                    # get satellite positions, SFR, and stellar masses
                    self.sub_pos = cat_sat.SubhaloPos
                    self.sub_SFR = gal_SFR
                    self.sub_sm = np.array(galprop['stellar_totmass'])



            for i in np.arange(self.n_selected_groups):
                if (i % size == rank):
                    print "Projecting i / task:", i, rank
                    grp_id = self.grp_ids[i]
                    gpos  = self.grp_pos[i]
                    r200  = self.grp_Rvir[i]
                    
                    if self.modify_temp: 
                        self.grid = np.zeros([self.n_cloudy_var,self.n_tempfac,self.n_species,self.ngrid,self.ngrid])
                    else:
                        self.grid = np.zeros([self.n_cloudy_var,self.n_species,self.ngrid,self.ngrid])
                    

                    CGMsnap_file_path = self.loadbase+ "s{}/{}.hdf5".format(snapnum,str(int(grp_id)).zfill(5))
                    outputpath = self.savebase+ "s{}/".format(snapnum)
                    savename = outputpath+str(int(grp_id)).zfill(5)+'.hdf5'

                    if not os.path.isfile(CGMsnap_file_path):
                        print "File {} does not exist!  Skipping...".format(CGMsnap_file_path)
                    elif os.path.isfile(savename):
                        print "File {} already exists!  Skipping...".format(savename)
                    else:
                        print "working on {}".format(savename)
                        [m,pos,metals,rho,u,nelec,hsml,T,neut_frac] = self.load_CGM_file(i,CGMsnap_file_path)

                        pos_cent = self._fix_pos(self.grp_pos[i],pos,self.box)
                        print "pos_cent: ",pos_cent
                        self.calc_grid(i,m,pos_cent,metals,rho,hsml,T,neut_frac,kernel_type='SPH')

                        print "Saving group file now to {}".format(savename)
                        print ""
                        # do not overwrite mode:
                        f=h5py.File(savename,'w-')
                        #over-write mode:
                        #f=h5py.File(outputpath+str(int(grp_id)).zfill(5)+'.hdf5','w')

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

                        if self.modify_temp:
                            for cc in xrange(self.n_cloudy_var):
                                cloudy = self.cloudy_list[cc]
                                c_grp = f.create_group(cloudy)
                                for tt in xrange(self.n_tempfac):
                                    tempfac = self.tempfac_arr[tt]
                                    t_grp = c_grp.create_group("temp_fac_{}".format(tempfac))
                                    for ss in xrange(self.n_species):
                                        elem = self.elem_list[ss]
                                        ion = self.ion_list[ss]
                                        dataset_name = elem+str(ion)
                                        t_grp.create_dataset(dataset_name,data=self.grid[cc,tt,ss,:,:])
                        elif not self.modify_temp:
                            for cc in xrange(self.n_cloudy_var):
                                cloudy = self.cloudy_list[cc]
                                c_grp = f.create_group(cloudy)
                                for ss in xrange(self.n_species):
                                    elem = self.elem_list[ss]
                                    ion = self.ion_list[ss]
                                    dataset_name = elem+str(ion)
                                    c_grp.create_dataset(dataset_name,data=self.grid[cc,ss,:,:])
                        f.close()




    def load_CGM_file(self,i,CGMsnap_file_path,use_block_name=True):
        f=h5py.File(CGMsnap_file_path,'r')
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


    def calc_grid(self,i,m,pos,metals,rho,hsml,T,neut_frac,kernel_type='SPH'):
        print "initial # of particles before picking out: ",np.size(m)
        for d in np.arange(3):
            ind = np.abs(pos[:,d]) <= 1.5*self.grid_radius
            pos = pos[ind]
            m = m[ind]
            metals = metals[ind]
            hsml = hsml[ind]
            T = T[ind]
            rho = rho[ind]  
            neut_frac = neut_frac[ind]
            #print "np.sum(ind) ",np.sum(ind)

        print "# of particles to fieldize over: ",np.size(m)
        r = AU.PhysicalPosition(np.sqrt(pos[:,0]**2.+pos[:,1]**2.+pos[:,2]**2.),self.redshift)

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
            r = r[cloudy_cut]

        pos_orig = np.copy(pos)
        m_orig = np.copy(m)
        metals_orig = np.copy(metals)
        hsml_orig = np.copy(hsml)
        T_orig = np.copy(T)
        rho_orig = np.copy(rho)
        rho_Hatoms_orig = np.copy(rho_Hatoms)
        neut_frac_orig = np.copy(neut_frac)
        r_orig = np.copy(r)

        for cc in xrange(self.n_cloudy_var):
            tab = self.tab_list[cc]
            cloudy = self.cloudy_list[cc]
            if self.modify_temp:
                for tt in xrange(self.n_tempfac):
                    tempfac = self.tempfac_arr[tt]
                    low_dens = np.copy(np.log10(rho_Hatoms) < -3.)
                    T = np.copy(T_orig)
                    T[low_dens] *= tempfac

                    T_cut = np.logical_and(np.log10(T) >= 3.,np.log10(T) <= 8.6)
                    pos = np.copy(pos_orig[T_cut])
                    m = np.copy(m_orig[T_cut])
                    metals = np.copy(metals_orig[T_cut])
                    hsml = np.copy(hsml_orig[T_cut])
                    T = np.copy(T[T_cut])
                    rho = np.copy(rho_orig[T_cut])
                    rho_Hatoms = np.copy(rho_Hatoms_orig[T_cut])
                    neut_frac = np.copy(neut_frac_orig[T_cut])
                    r = np.copy(r_orig[T_cut])

                    for ss in xrange(self.n_species):
                        elem = self.elem_list[ss]
                        ion = self.ion_list[ss]
                        print "species: {}".format(elem+str(ion))
                        print "cloudy: {}".format(cloudy)
                        print "tempfac: {}".format(tempfac)

                        if elem == "H" and ion == 1:
                            H_massfrac = metals[:,0]
                            star=cold_gas.RahmatiRT(self.redshift, self.hubble)
                            fake_bar = {'Density':rho,'NeutralHydrogenAbundance':neut_frac}
                            new_neut_frac = star.get_reproc_HI(fake_bar)
                            species_mass = m*H_massfrac * new_neut_frac

                        else:
                            elem_massfrac = metals[:,AU.elem_lookup(elem)]
                            if cloudy == "ion_out_fancy_atten":
                                ion_frac = tab.ion(elem,ion,rho_Hatoms,T)
                            elif cloudy == "UVB_sf_xrays_ext":
                                beta = self.grp_SFR[i]/r**2.
                                print "SFR = ",self.grp_SFR[i]
                                ion_frac = tab.ion(elem,ion,rho_Hatoms,T,beta=beta)
                            species_mass = m*elem_massfrac*ion_frac

                        temp_grid = np.zeros([self.ngrid,self.ngrid])
                        self.sub_gridize_single_halo(pos,hsml,species_mass,temp_grid,kernel_type=kernel_type)
                        self.grid[cc,tt,ss,:,:] = temp_grid

                        elem_mass_g = AU.elem_atom_mass(elem)
                        massg=AU.UnitMass_in_g/self.hubble*(1/elem_mass_g)
                        epsilon=2.*self.grid_radius/self.ngrid*AU.UnitLength_in_cm/self.hubble/(1+self.redshift)
                        self.grid[cc,tt,ss,:,:] *= (massg/epsilon**2)
                        self.grid[cc,tt,ss,:,:] += 0.1
                        np.log10(self.grid[cc,tt,ss,:,:],self.grid[cc,tt,ss,:,:])


            elif not self.modify_temp:
                for ss in xrange(self.n_species):
                    elem = self.elem_list[ss]
                    ion = self.ion_list[ss]
                    print "species: {}".format(elem+str(ion))
                    print "cloudy: {}".format(cloudy)

                    if elem == "H" and ion == 1:
                        H_massfrac = metals[:,0]
                        star=cold_gas.RahmatiRT(self.redshift, self.hubble)
                        fake_bar = {'Density':rho,'NeutralHydrogenAbundance':neut_frac}
                        new_neut_frac = star.get_reproc_HI(fake_bar)
                        species_mass = m*H_massfrac * new_neut_frac

                    else:
                        elem_massfrac = metals[:,AU.elem_lookup(elem)]
                        if cloudy == "ion_out_fancy_atten":
                            ion_frac = tab.ion(elem,ion,rho_Hatoms,T)
                        elif cloudy == "UVB_sf_xrays_ext":
                            if not self.multiple_sources:
                                beta = self.grp_SFR[i]/r**2.
                                print "SFR = ",self.grp_SFR[i]
                            elif self.multiple_sources:
                                grp_id = self.grp_ids[i]
                                sub_ids = self.subhalo_id_dict[grp_id]
                                n_sources = np.size(sub_ids)
                                if n_sources == 0: 
                                    print "No subhalos for group {}!  ERROR".format(grp_id)
                                elif n_sources == 1: #only central galaxy
                                    beta = self.grp_SFR[i]/r**2.
                                else:
                                    # loop over all sources; for each source, get its position, SFR, and stellar mass
                                    beta_arr = np.zeros([n_sources,np.size(r)])
                                    for kk in np.arange(n_sources):
                                        sub_id = sub_ids[kk]
                                        sub_SFR = self.sub_SFR[sub_id]
                                        sub_sm = self.sub_sm[sub_id]

                                        # HERE: only grid if sm is above some min threshold mass
                                        if sub_SFR == 0.:
                                            pass # leave beta = 0 for this source everywhere
                                        else:
                                            # "Fix position" of source, same way that we did with particles:
                                            subpos_cent = self._fix_pos(self.grp_pos[i],np.array([self.sub_pos[sub_id]]),self.box)
                                            subpos_cent = subpos_cent[0]
                                            # Then calculate distance from each cell to the subhalo position
                                            rs = AU.PhysicalPosition(np.sqrt(\
                                            (pos[:,0]-subpos_cent[0])**2.+\
                                            (pos[:,1]-subpos_cent[1])**2.+\
                                            (pos[:,2]-subpos_cent[2])**2.)\
                                            ,self.redshift)
                                            
                                            if np.sum(rs > 2000.) > 0: print "radius error!"

                                            beta_arr[kk,:] = sub_SFR/rs**2.
                                    beta = np.sum(beta_arr,axis=0)
                                    print "beta ",beta
                                    print "np.shape(beta) ",np.shape(beta)
                                    print "np.shape(r) ",np.shape(r)

                            # Given beta for each cell (radiation that each cell sees from young stellar pops), 
                            # calculate ionization states
                            ion_frac = tab.ion(elem,ion,rho_Hatoms,T,beta=beta)

                        species_mass = m*elem_massfrac*ion_frac

                    temp_grid = np.zeros([self.ngrid,self.ngrid])
                    self.sub_gridize_single_halo(pos,hsml,species_mass,temp_grid,kernel_type=kernel_type)
                    self.grid[cc,ss,:,:] = temp_grid

                    elem_mass_g = AU.elem_atom_mass(elem)
                    massg=AU.UnitMass_in_g/self.hubble*(1/elem_mass_g)
                    epsilon=2.*self.grid_radius/self.ngrid*AU.UnitLength_in_cm/self.hubble/(1+self.redshift)
                    self.grid[cc,ss,:,:] *= (massg/epsilon**2)
                    self.grid[cc,ss,:,:] += 0.1
                    np.log10(self.grid[cc,ss,:,:],self.grid[cc,ss,:,:])




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








