# Unlike CGM_extract, which extracts all gas within a specific radius of each desired halo, this
# script focuses on extracting gas that is bound to the FoF groups, thus exploiting the natural
# format of the snapshots themselves.

import numpy as np

import readsubfHDF5 
import snapHDF5 
import numpy as np
import os
import sys
import h5py
import time

from units import AREPO_units
AU = AREPO_units()

from mpi4py import MPI

import convert_cloudy as cc
import cold_gas


class illustris_CGM_extract:
    def __init__(self,snapnum):
        run = "ill1"
        self.snapnum = np.int32(snapnum)
        if self.snapnum == 135: self.redshift = 0.
        elif self.snapnum == 127: self.redshift = 0.1
        elif self.snapnum == 120: self.redshift = 0.2
        elif self.snapnum == 103: self.redshift = 0.5
        elif self.snapnum == 85: self.redshift = 1.
        elif self.snapnum == 68: self.redshift = 2.
        elif self.snapnum == 60: self.redshift = 3.
        elif self.snapnum == 54: self.redshift = 4.

        print "self.redshift ",self.redshift

        self.hubble = 0.7
        # self.snapnum_arr = np.array([135,68,127,120,103,85,60,54])])
        # self.redshift_arr = False #NOT IMPLEMENTED 

        verbose=False

        group_min_mass = 10.**11
        group_max_mass = 10.**16
        dat_str_list = ["Masses","GFM_Metals","GFM_Metallicity","Density","Volume","InternalEnergy","ElectronAbundance","NeutralHydrogenAbundance"]
        self.block_list = ["MASS","GMET","GZ  ","RHO ","VOL ","U   ","NE  ","NH  "]
        self.savebase = '/n/home04/jsuresh/scratch1/CGM_prop/'

        if run == 'ill3': self.snapbase="/n/ghernquist/Illustris/Runs/Illustris-3/output/"
        elif run == 'ill2': self.snapbase="/n/ghernquist/Illustris/Runs/Illustris-2/output/"
        elif run == 'ill1': self.snapbase="/n/ghernquist/Illustris/Runs/Illustris-1/output/"
        # self.snapbase = "/n/hernquistfs1/Illustris/SmallBox/GFM/Production/Cosmo/Cosmo0_V6/L25n256/output/"

        self.tab = cc.CloudyTable(self.redshift,atten_type="ion_out_fancy_atten")

        # Here's where the magic happens
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        print "my rank = {}".format(rank)
        size = comm.Get_size()
        print "my size = {}".format(size)

        if rank == 0:
            [self.Group_M_Crit200,self.SubhaloLenType,self.GroupLenType,self.GroupFirstSub,\
            self.GroupOffset,self.HaloOffset,self.FileNum,self.FileTypeNumbers,self.multiple]\
            = self.setup_tables(self.snapnum,verbose=verbose)
        
            self.gm = AU.PhysicalMass(self.Group_M_Crit200)
            self.grnr_arr = np.arange(np.size(self.gm))
            mass_cut = np.logical_and(self.gm > group_min_mass,self.gm < group_max_mass)
            self.grnr_arr = self.grnr_arr[mass_cut]
            self.gm = self.gm[mass_cut]
        else: 
            self.Group_M_Crit200 = None
            self.SubhaloLenType = None
            self.GroupLenType = None
            self.GroupFirstSub = None
            self.GroupOffset = None
            self.HaloOffset = None
            self.FileNum = None
            self.FileTypeNumbers = None
            self.multiple = None
            self.grnr_arr = None

        # self.Group_M_Crit200 = comm.bcast(self.Group_M_Crit200,root=0)
        self.SubhaloLenType = comm.bcast(self.SubhaloLenType,root=0)
        self.GroupLenType = comm.bcast(self.GroupLenType,root=0)
        self.GroupFirstSub = comm.bcast(self.GroupFirstSub,root=0)
        self.GroupOffset = comm.bcast(self.GroupOffset,root=0)
        self.HaloOffset = comm.bcast(self.HaloOffset,root=0)
        self.FileNum = comm.bcast(self.FileNum,root=0)
        self.FileTypeNumbers = comm.bcast(self.FileTypeNumbers,root=0)
        self.multiple = comm.bcast(self.multiple,root=0)
        self.grnr_arr = comm.bcast(self.grnr_arr,root=0)

        
        n_grp = np.size(self.grnr_arr)

        if rank == 0:
            self.phase_mass = np.zeros([n_grp,10]) #include ISM and CGM as first two phases
            self.phase_metallicity = np.zeros([n_grp,10])
            self.phase_totalMetalMass = np.zeros([n_grp,10]) 
            self.phase_metalsMass = np.zeros([n_grp,10,9]) 
            self.phase_volumeFrac = np.zeros([n_grp,10])
            self.phase_HI_mass = np.zeros([n_grp,10])
            self.phase_MgII_mass = np.zeros([n_grp,10])
            self.phase_CIII_mass = np.zeros([n_grp,10])
            self.phase_CIV_mass = np.zeros([n_grp,10])
            self.phase_SiIII_mass = np.zeros([n_grp,10])
            self.phase_SiIV_mass = np.zeros([n_grp,10])
            self.phase_OVI_mass = np.zeros([n_grp,10])
            self.phase_OVII_mass = np.zeros([n_grp,10])
            self.phase_OVIII_mass = np.zeros([n_grp,10])
            self.phase_NeVIII_mass = np.zeros([n_grp,10])
            self.phase_FeII_mass = np.zeros([n_grp,10])


        for i in np.arange(n_grp):
            rank_to_work_on_this_grp = (i % (size-1))+1
            if (rank == rank_to_work_on_this_grp):
                print "Projecting i / task:", i, rank
                print "On group {} of {}".format(i+1,n_grp)
                ts = time.time()
                grnr = self.grnr_arr[i]
                data_dict = self.read_data(grnr,verbose=verbose)
                CGM_dict = self.calc_desired_quantities(grnr,data_dict)
                print "Time for group {}: {}".format(grnr,time.time()-ts)

                print "sending phase_mass: ",CGM_dict['phase_mass']
                comm.send(obj=CGM_dict,dest=0,tag=i)
                del CGM_dict

            if rank == 0:
                CGM_dict = comm.recv(source=rank_to_work_on_this_grp,tag=i)
                self.phase_mass[i,:] = CGM_dict['phase_mass']
                self.phase_metallicity[i,:] = CGM_dict['phase_metallicity']
                self.phase_totalMetalMass[i,:] = CGM_dict['phase_totalMetalMass']
                self.phase_metalsMass[i,:] = CGM_dict['phase_metalsMass']
                self.phase_volumeFrac[i,:] = CGM_dict['phase_volumeFrac']
                self.phase_HI_mass[i,:] = CGM_dict['phase_HI_mass']
                self.phase_MgII_mass[i,:] = CGM_dict['phase_MgII_mass']
                self.phase_CIII_mass[i,:] = CGM_dict['phase_CIII_mass']
                self.phase_CIV_mass[i,:] = CGM_dict['phase_CIV_mass']
                self.phase_SiIII_mass[i,:] = CGM_dict['phase_SiIII_mass']
                self.phase_SiIV_mass[i,:] = CGM_dict['phase_SiIV_mass']
                self.phase_OVI_mass[i,:] = CGM_dict['phase_OVI_mass']
                self.phase_OVII_mass[i,:] = CGM_dict['phase_OVII_mass']
                self.phase_NeVIII_mass[i,:] = CGM_dict['phase_NeVIII_mass']
                self.phase_FeII_mass[i,:] = CGM_dict['phase_FeII_mass']


        if rank == 0:
            print "self.phase_mass ",self.phase_mass
            print "self.phase_mass[:,0] ",self.phase_mass[:,0]
            self.save_data(verbose=verbose)
                



    def setup_tables(self,snapnum,parttype=0,long_ids=True,double_output=False,verbose=False):
        # global FlagRead, cat, GroupOffset, HaloOffset, multiple, filename, Parttype, FileTypeNumbers, FileNum
 
        print "Setting up tables..."
        ts = time.time()

        if (verbose):
            print "READHALO: INITIAL READ"
            print "READHALO: Parttype = ", parttype

        #read in catalog
        cat = readsubfHDF5.subfind_catalog(self.snapbase, self.snapnum, long_ids=long_ids, double_output=double_output, keysel=["GroupLenType","GroupNsubs","GroupFirstSub","SubhaloLenType","SubhaloMassType","SubhaloGrNr","Group_M_Crit200","GroupPos","Group_R_Crit200"])
        Group_M_Crit200 = np.copy(cat.Group_M_Crit200)
        SubhaloLenType = np.copy(cat.SubhaloLenType)
        GroupLenType = np.copy(cat.GroupLenType)
        GroupFirstSub = np.copy(cat.GroupFirstSub)

        if (cat.ngroups==0):
            if (verbose):
                print "READHALO: no groups in catalog... returning"
            return

        #Construct offset tables
        GroupOffset = np.zeros([cat.ngroups, 6], dtype="int64")
        HaloOffset = np.zeros([cat.nsubs, 6], dtype="int64")
        filename = self.snapbase+"/snap_"+str(self.snapnum).zfill(3)
        multiple=False
        if (os.path.exists(filename+".hdf5")==False):
                filename = self.snapbase+"/snapdir_"+str(self.snapnum).zfill(3)+"/snap_"+str(self.snapnum).zfill(3)+"."+str(0)
                multiple=True
        if (os.path.exists(filename+".hdf5")==False):
                print "READHALO: [error] file not found : ", filename
                sys.exit()

        k=0
        for i in range(0, cat.ngroups):
            if (i>0):
                GroupOffset[i, parttype] =  GroupOffset[i-1, parttype] + cat.GroupLenType[i-1, parttype]
            if (cat.GroupNsubs[i]>0):
                HaloOffset[k, parttype] = GroupOffset[i, parttype]
                k+=1
                for j in range(1, cat.GroupNsubs[i]):
                    HaloOffset[k, parttype] =  HaloOffset[k-1, parttype] + cat.SubhaloLenType[k-1, parttype]
                    k+=1
        if (k!=cat.nsubs):
            print "READHALO: problem with offset table", k, cat.nsubs
            sys.exit()

        #construct file tables
        if (multiple):
            filename = self.snapbase+"/snapdir_"+str(self.snapnum).zfill(3)+"/snap_"+str(self.snapnum).zfill(3)+"."+str(0)
        else:
            filename = self.snapbase+"/snap_"+str(self.snapnum).zfill(3)

        head = snapHDF5.snapshot_header(filename)
        FileNum = head.filenum

        FileTypeNumbers = np.zeros([FileNum, 6], dtype="int64") 
        cumcount = np.zeros(6, dtype="int64")

        for fnr in range(0, FileNum-1):
            if (multiple):
                filename = self.snapbase+"/snapdir_"+str(self.snapnum).zfill(3)+"/snap_"+str(self.snapnum).zfill(3)+"."+str(fnr)
            else:
                filename = self.snapbase+"/snap_"+str(self.snapnum).zfill(3)

            if (verbose):
                print "READHALO: initial reading file :", filename

            head = snapHDF5.snapshot_header(filename)
    
            cumcount[:] += head.npart[:]
            FileTypeNumbers[fnr+1, :] = cumcount[:]

        print "Time to set up tables: ",time.time()-ts

        return [Group_M_Crit200,SubhaloLenType,GroupLenType,GroupFirstSub,GroupOffset,HaloOffset,FileNum,FileTypeNumbers,multiple]



    def read_data(self,fof_num,sub_num=-1,parttype=0,verbose=False):
        data_dict = {}

        if (sub_num>=0) & (fof_num < 0):
            off = HaloOffset[sub_num, parttype]
            left = SubhaloLenType[sub_num, parttype]
            # if (verbose):
            #     print "READHALO: nr / particle # / mass :", sub_num, self.cat.SubhaloLenType[sub_num, parttype], self.cat.SubhaloMassType[sub_num, parttype].astype("float64")
        if (fof_num>=0) & (sub_num < 0):
            off = self.GroupOffset[fof_num, parttype]
            left = self.GroupLenType[fof_num, parttype]
            # if (verbose):
            #     print "READHALO: nr / particle # / mass :", fof_num, self.cat.GroupLenType[fof_num, parttype], self.cat.GroupMassType[fof_num, parttype].astype("float64")
        if (sub_num>=0) & (fof_num>=0):
            real_sub_num = sub_num + GroupFirstSub[fof_num]
            off = HaloOffset[real_sub_num, parttype]
            left = SubhaloLenType[real_sub_num, parttype]
            # if (verbose):
            #     print "READHALO: nr / particle # / mass :", real_sub_num, SubhaloLenType[real_sub_num, parttype], SubhaloMassType[real_sub_num, parttype].astype("float64")
        
        if (left==0 and verbose):
            print "READHALO: no particles of type... returning"
            return

        # Get first file that contains particles of required halo/fof/etc
        findex = np.argmax(self.FileTypeNumbers[:, parttype] > off) - 1
        #in case we reached the end argmax returns 0
        if (findex == -1): 
            findex = self.FileNum - 1

        if (verbose):
            print "READHALO: first file that contains particles =", findex

        for fnr in range(0, findex):
            off -= self.FileTypeNumbers[fnr+1, parttype] - self.FileTypeNumbers[fnr, parttype] 


        # Read data from file
        first=True
        for fnr in range(findex, self.FileNum):
            if (self.multiple):
                filename = self.snapbase+"/snapdir_"+str(self.snapnum).zfill(3)+"/snap_"+str(self.snapnum).zfill(3)+"."+str(fnr)
            else:
                filename = self.snapbase+"/snap_"+str(self.snapnum).zfill(3)
        
            if (verbose):
                print "READHALO: reading file :", filename
        
            head = snapHDF5.snapshot_header(filename)
            nloc = head.npart[parttype]

            if (nloc > off):
                if (verbose):
                    print "READHALO: data"
                start = off
                if (nloc - off > left):
                    count = left    
                else:
                    count = nloc - off

                if (first==True):   
                    for block_name in self.block_list:
                        if (verbose):
                            print "Reading block: ",block_name
                        data_dict[block_name] = snapHDF5.read_block(filename, block_name, parttype, slab_start=start, slab_len=count)
                    first=False
                else:
                    for block_name in self.block_list:
                        if (verbose): print "Reading block again: ",block_name
                        # print "Reading block again: ",block_name
                        print "start ",start
                        print "len ",count
                        data_dict[block_name] = np.append(data_dict[block_name], snapHDF5.read_block(filename, block_name, parttype, slab_start=start, slab_len=count), axis=0)
                left -= count
                off += count
            if (left==0):
                break
            off -= nloc

        return data_dict



    def calc_desired_quantities(self,grnr,data_dict):
        # self.block_list = ["MASS","GMET","GZ  ","RHO ","VOL ","U   ","NE  ","NH  "]
        # Determine which gas is in each phase:
        rho = AU.PhysicalDensity(data_dict["RHO "],self.redshift)
        u = np.array(data_dict["U   "])
        nelec = np.array(data_dict["NE  "])
        T = AU.GetTemp(u, nelec, gamma = 5.0/3.0)
        z = data_dict["GZ  "]
        met = data_dict["GMET"]
        vol = data_dict["VOL "]
        m = AU.PhysicalMass(data_dict["MASS"])
        v_tot = np.sum(vol)

        # HI:
        H_massfrac = met[:,0]
        star=cold_gas.RahmatiRT(self.redshift, self.hubble)
        fake_bar = {'Density':data_dict["RHO "],'NeutralHydrogenAbundance':data_dict["NH  "]}
        new_neut_frac = star.get_reproc_HI(fake_bar)
        HI_mass = m*H_massfrac * new_neut_frac
        # Other ions:
        # Prep cloudy
        T_cut = np.logical_and(np.log10(T) >= 3.,np.log10(T) <= 8.6)
        rho_Hatoms = rho*(met[:,0]/AU.ProtonMass)
        rho_cut = np.logical_and(np.log10(rho_Hatoms) >= -7.,np.log10(rho_Hatoms) <= 4.)
        cloudy_cut = np.logical_and(T_cut,rho_cut)
        print "Cloudy cut removed {} particles ".format(np.sum(np.logical_not(cloudy_cut)))
        rho = rho[cloudy_cut]
        u = u[cloudy_cut]
        nelec = nelec[cloudy_cut]
        T = T[cloudy_cut]
        z = z[cloudy_cut]
        met = met[cloudy_cut]
        vol = vol[cloudy_cut]
        m = m[cloudy_cut]
        v_tot = np.sum(vol)
        rho_Hatoms = rho_Hatoms[cloudy_cut]


        HI_mass = np.zeros(0)
        MgII_mass = np.zeros(0)
        CIII_mass = np.zeros(0)
        CIV_mass = np.zeros(0)
        SiIII_mass = np.zeros(0)
        SiIV_mass = np.zeros(0)
        OI_mass = np.zeros(0)
        OVI_mass = np.zeros(0)
        OVII_mass = np.zeros(0)
        OVIII_mass = np.zeros(0)
        NeVIII_mass = np.zeros(0)
        FeII_mass = np.zeros(0)
        elem_list = ['H','Mg','C','C','Si','Si','O','O','O','Ne','Fe']
        ion_list = [1,2,3,4,3,4,6,7,8,8,2]
        # arr_list = [HI_mass,MgII_mass,CIII_mass,CIV_mass,SiIII_mass,SiIV_mass,OVI_mass,OVII_mass,OVIII_mass,NeVIII_mass,FeII_mass]
        n_spec = len(elem_list)
        for k in xrange(n_spec):
            elem = elem_list[k]
            ion = ion_list[k]
            elem_massfrac = met[:,AU.elem_lookup(elem)]
            ion_frac = self.tab.ion(elem,ion,rho_Hatoms,T)
            species_mass = m*elem_massfrac*ion_frac
            if k == 0: HI_mass = species_mass
            elif k == 1: MgII_mass = species_mass
            elif k == 2: CIII_mass = species_mass
            elif k == 3: CIV_mass = species_mass
            elif k == 4: SiIII_mass = species_mass
            elif k == 5: SiIV_mass = species_mass
            elif k == 6: OVI_mass = species_mass
            elif k == 7: OVII_mass = species_mass
            elif k == 8: OVIII_mass = species_mass
            elif k == 9: NeVIII_mass = species_mass
            elif k == 10: FeII_mass = species_mass



        m_typ = 2.1e-24 # mass of typical particle (75% H, 25% He)
        rho /= m_typ

        ISM = rho > 0.13
        CGM = np.logical_not(ISM)
        cold_den = np.logical_and(CGM,np.logical_and(rho>10.**-3.,T<=10.**4.))
        cool_den = np.logical_and(CGM,np.logical_and(rho>10.**-3.,np.logical_and(T>10.**4.,T<=10.**5.)))
        warm_den = np.logical_and(CGM,np.logical_and(rho>10.**-3.,np.logical_and(T>10.**5.,T<=10.**6.)))
        hot_den = np.logical_and(CGM,np.logical_and(T>10.**6.,rho>10.**-3.))
        cold_dif = np.logical_and(CGM,np.logical_and(rho<=10.**-3.,T<=10.**4.))
        cool_dif = np.logical_and(CGM,np.logical_and(rho<=10.**-3.,np.logical_and(T>10.**4.,T<=10.**5.)))
        warm_dif = np.logical_and(CGM,np.logical_and(rho<=10.**-3.,np.logical_and(T>10.**5.,T<=10.**6.)))
        hot_dif = np.logical_and(CGM,np.logical_and(T>10.**6.,rho<=10.**-3.))

        phase_bool_list = [ISM,CGM,cold_den,cool_den,warm_den,hot_den,cold_dif,cool_dif,warm_dif,hot_dif]
        n_phase = len(phase_bool_list)

        # for j in xrange(n_phase):
        #     phase = phase_bool_list[j]
        #     self.phase_mass[i,j] = np.sum(m[phase])
        #     self.phase_metallicity[i,j] = np.sum(m[phase]*z[phase])/self.phase_mass[i,j]
        #     self.phase_totalMetalMass[i,j] = self.phase_mass[i,j]*self.phase_metallicity[i,j]
        #     for k in np.arange(9):
        #         self.phase_metalsMass[i,j,k] = np.sum(m[phase]*met[phase,k])
        #     self.phase_volumeFrac[i,j] = np.sum(vol[phase])/v_tot
        #     self.phase_HI_mass[i,j] = np.sum(HI_mass[phase])
        #     self.phase_MgII_mass[i,j] = np.sum(MgII_mass[phase])
        #     self.phase_CIII_mass[i,j] = np.sum(CIII_mass[phase])
        #     self.phase_CIV_mass[i,j] = np.sum(CIV_mass[phase])
        #     self.phase_SiIII_mass[i,j] = np.sum(SiIII_mass[phase])
        #     self.phase_SiIV_mass[i,j] = np.sum(SiIV_mass[phase])
        #     self.phase_OVI_mass[i,j] = np.sum(OVI_mass[phase])
        #     self.phase_OVII_mass[i,j] = np.sum(OVII_mass[phase])
        #     self.phase_OVIII_mass[i,j] = np.sum(OVIII_mass[phase])
        #     self.phase_NeVIII_mass[i,j] = np.sum(NeVIII_mass[phase])
        #     self.phase_FeII_mass[i,j] = np.sum(FeII_mass[phase])


        phase_mass = np.zeros(10)
        phase_metallicity = np.zeros(10)
        phase_totalMetalMass = np.zeros(10) 
        phase_metalsMass = np.zeros([10,9]) 
        phase_volumeFrac = np.zeros(10)
        phase_HI_mass = np.zeros(10)
        phase_MgII_mass = np.zeros(10)
        phase_CIII_mass = np.zeros(10)
        phase_CIV_mass = np.zeros(10)
        phase_SiIII_mass = np.zeros(10)
        phase_SiIV_mass = np.zeros(10)
        phase_OVI_mass = np.zeros(10)
        phase_OVII_mass = np.zeros(10)
        phase_OVIII_mass = np.zeros(10)
        phase_NeVIII_mass = np.zeros(10)
        phase_FeII_mass = np.zeros(10)

        for j in xrange(n_phase):
            phase = phase_bool_list[j]
            phase_mass[j] = np.sum(m[phase])
            phase_metallicity[j] = np.sum(m[phase]*z[phase])/phase_mass[j]
            phase_totalMetalMass[j] = phase_mass[j]*phase_metallicity[j]
            for k in np.arange(9):
                phase_metalsMass[j,k] = np.sum(m[phase]*met[phase,k])
            phase_volumeFrac[j] = np.sum(vol[phase])/v_tot
            phase_HI_mass[j] = np.sum(HI_mass[phase])
            phase_MgII_mass[j] = np.sum(MgII_mass[phase])
            phase_CIII_mass[j] = np.sum(CIII_mass[phase])
            phase_CIV_mass[j] = np.sum(CIV_mass[phase])
            phase_SiIII_mass[j] = np.sum(SiIII_mass[phase])
            phase_SiIV_mass[j] = np.sum(SiIV_mass[phase])
            phase_OVI_mass[j] = np.sum(OVI_mass[phase])
            phase_OVII_mass[j] = np.sum(OVII_mass[phase])
            phase_OVIII_mass[j] = np.sum(OVIII_mass[phase])
            phase_NeVIII_mass[j] = np.sum(NeVIII_mass[phase])
            phase_FeII_mass[j] = np.sum(FeII_mass[phase])

        CGM_dict = {'grnr':grnr,'phase':phase,'phase_mass':phase_mass,'phase_metallicity':phase_metallicity,\
        'phase_totalMetalMass':phase_totalMetalMass,'phase_metalsMass':phase_metalsMass,'phase_volumeFrac':phase_volumeFrac,\
        'phase_HI_mass':phase_HI_mass,'phase_MgII_mass':phase_MgII_mass,'phase_CIII_mass':phase_CIII_mass,\
        'phase_CIV_mass':phase_CIV_mass,'phase_SiIII_mass':phase_SiIII_mass,'phase_SiIV_mass':phase_SiIV_mass,\
        'phase_OVI_mass':phase_OVI_mass,'phase_OVII_mass':phase_OVII_mass,'phase_OVIII_mass':phase_OVIII_mass,\
        'phase_NeVIII_mass':phase_NeVIII_mass,'phase_FeII_mass':phase_FeII_mass}

        return CGM_dict


    def save_data(self,verbose=False):
        filepath = self.savebase+"s{}.hdf5".format(self.snapnum)

        if verbose:
            print "Saving CGM snapshot for snapshot {} to {}".format(self.snapnum,filepath)

        f=h5py.File(filepath,'a')

        head = f.create_group("Header")
        head.create_dataset("grp_id",data=self.grnr_arr)
        head.create_dataset("grp_mass",data=self.gm)

        p_grp = f.create_group("Gas")
        p_grp.create_dataset("phase_mass",data=self.phase_mass)
        p_grp.create_dataset("phase_metallicity",data=self.phase_metallicity)
        p_grp.create_dataset("phase_totalMetalMass",data=self.phase_totalMetalMass)
        p_grp.create_dataset("phase_metalsMass",data=self.phase_metalsMass)
        p_grp.create_dataset("phase_volumeFrac",data=self.phase_volumeFrac)
        p_grp.create_dataset("phase_HI_mass",data=self.phase_HI_mass)
        p_grp.create_dataset("phase_MgII_mass",data=self.phase_MgII_mass)
        p_grp.create_dataset("phase_CIII_mass",data=self.phase_CIII_mass)
        p_grp.create_dataset("phase_CIV_mass",data=self.phase_CIV_mass)
        p_grp.create_dataset("phase_SiIII_mass",data=self.phase_SiIII_mass)
        p_grp.create_dataset("phase_SiIV_mass",data=self.phase_SiIV_mass)
        p_grp.create_dataset("phase_OVI_mass",data=self.phase_OVI_mass)
        p_grp.create_dataset("phase_OVII_mass",data=self.phase_OVII_mass)
        p_grp.create_dataset("phase_OVIII_mass",data=self.phase_OVIII_mass)
        p_grp.create_dataset("phase_NeVIII_mass",data=self.phase_NeVIII_mass)
        p_grp.create_dataset("phase_FeII_mass",data=self.phase_FeII_mass)

        f.close()



if __name__ == '__main__':
    for arg in sys.argv:
        print arg
    illustris_CGM_extract(arg)
