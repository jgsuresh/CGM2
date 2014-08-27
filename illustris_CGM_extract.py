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



class illustris_CGM_extract:
    def __init__(self,verbose=False):
        run = "ill1"
        self.snapnum = 120

        # group_min_mass = 10.**11.8 
        # group_max_mass = 10.**12.2
        gal_min_stellar_mass = 10.**10.
        gal_max_stellar_mass = 10.**11.6
        dat_str_list = ["Masses","Coordinates","GFM_Metals","GFM_Metallicity","Velocities","Density","Volume","InternalEnergy","ElectronAbundance","NeutralHydrogenAbundance","SmoothingLength"]
        self.block_list = ["MASS","POS ","GMET","GZ  ","VEL ","RHO ","VOL ","U   ","NE  ","NH  ","HSML"]
        self.savebase = '/n/home04/jsuresh/scratch1/AREPOfest/data/CGM_snaps/'

        if run == 'ill3': self.snapbase="/n/ghernquist/Illustris/Runs/Illustris-3/output/"
        elif run == 'ill2': self.snapbase="/n/ghernquist/Illustris/Runs/Illustris-2/output/"
        elif run == 'ill1': self.snapbase="/n/ghernquist/Illustris/Runs/Illustris-1/output/"

        self.setup_tables(verbose=verbose)

        grnr_arr = self.get_grnr_array(gal_min_stellar_mass,gal_max_stellar_mass)
        # grnr_arr = np.array([79])

        for grnr in grnr_arr:
            ts = time.time()
            filepath = self.savebase+"s{}/{}.hdf5".format(self.snapnum,str(int(grnr)).zfill(5))
            if os.path.isfile(filepath):
                print "File {} already exists!  Skipping...".format(filepath)
            else:
                data_dict = self.read_data(grnr,verbose=verbose)
                self.save_data(grnr,data_dict,verbose=verbose)
                print "Time for group {}: {}".format(grnr,time.time()-ts)




    def setup_tables(self,parttype=0,long_ids=True,double_output=False,verbose=False):
        # global FlagRead, cat, GroupOffset, HaloOffset, multiple, filename, Parttype, FileTypeNumbers, FileNum
 
        print "Setting up tables..."
        ts = time.time()

        if (verbose):
            print "READHALO: INITIAL READ"
            print "READHALO: Parttype = ", parttype

        #read in catalog
        self.cat = readsubfHDF5.subfind_catalog(self.snapbase, self.snapnum, long_ids=long_ids, double_output=double_output, keysel=["GroupLenType","GroupNsubs","GroupFirstSub","SubhaloLenType","SubhaloMassType","SubhaloGrNr","Group_M_Crit200","GroupPos","Group_R_Crit200"])

        if (self.cat.ngroups==0):
            if (verbose):
                print "READHALO: no groups in catalog... returning"
            return

        #Construct offset tables
        self.GroupOffset = np.zeros([self.cat.ngroups, 6], dtype="int64")
        self.HaloOffset = np.zeros([self.cat.nsubs, 6], dtype="int64")
        filename = self.snapbase+"/snap_"+str(self.snapnum).zfill(3)
        self.multiple=False
        if (os.path.exists(filename+".hdf5")==False):
                filename = self.snapbase+"/snapdir_"+str(self.snapnum).zfill(3)+"/snap_"+str(self.snapnum).zfill(3)+"."+str(0)
                self.multiple=True
        if (os.path.exists(filename+".hdf5")==False):
                print "READHALO: [error] file not found : ", filename
                sys.exit()

        k=0
        for i in range(0, self.cat.ngroups):
            if (i>0):
                self.GroupOffset[i, parttype] =  self.GroupOffset[i-1, parttype] + self.cat.GroupLenType[i-1, parttype]
            if (self.cat.GroupNsubs[i]>0):
                self.HaloOffset[k, parttype] = self.GroupOffset[i, parttype]
                k+=1
                for j in range(1, self.cat.GroupNsubs[i]):
                    self.HaloOffset[k, parttype] =  self.HaloOffset[k-1, parttype] + self.cat.SubhaloLenType[k-1, parttype]
                    k+=1
        if (k!=self.cat.nsubs):
            print "READHALO: problem with offset table", k, self.cat.nsubs
            sys.exit()

        #construct file tables
        if (self.multiple):
            filename = self.snapbase+"/snapdir_"+str(self.snapnum).zfill(3)+"/snap_"+str(self.snapnum).zfill(3)+"."+str(0)
        else:
            filename = self.snapbase+"/snap_"+str(self.snapnum).zfill(3)

        head = snapHDF5.snapshot_header(filename)
        self.FileNum = head.filenum

        self.FileTypeNumbers = np.zeros([self.FileNum, 6], dtype="int64") 
        cumcount = np.zeros(6, dtype="int64")

        for fnr in range(0, self.FileNum-1):
            if (self.multiple):
                filename = self.snapbase+"/snapdir_"+str(self.snapnum).zfill(3)+"/snap_"+str(self.snapnum).zfill(3)+"."+str(fnr)
            else:
                filename = self.snapbase+"/snap_"+str(self.snapnum).zfill(3)

            if (verbose):
                print "READHALO: initial reading file :", filename

            head = snapHDF5.snapshot_header(filename)
    
            cumcount[:] += head.npart[:]
            self.FileTypeNumbers[fnr+1, :] = cumcount[:]

        print "Time to set up tables: ",time.time()-ts



    def get_grnr_array(self,gal_min_stellar_mass,gal_max_stellar_mass):
        # Load galaxy properties catalog
        galf = h5py.File(self.snapbase+"/postprocessing/galprop/galprop_{}.hdf5".format(self.snapnum),'r')
        gal_sm = AU.PhysicalMass(np.array(galf['stellar_totmass']))
        n_gal = np.size(gal_sm)

        sm_mass_select = np.logical_and(gal_sm > gal_min_stellar_mass,gal_sm < gal_max_stellar_mass)
        gal_ids = np.arange(n_gal)[sm_mass_select]

        # filter out only primary halos:
        # method 1 - must be most massive stellar-wise in its group [not implemented yet]
        # method 2 - must inhabit most massive subhalo in the group
        primary_gal_ids = np.ones(0)
        for gal_id in gal_ids:
            grnr = self.cat.SubhaloGrNr[gal_id]
            if self.cat.GroupFirstSub[grnr] ==  gal_id: 
                primary_gal_ids = np.append(primary_gal_ids,gal_id)
        # print "np.size(primary_gal_ids) ",np.size(primary_gal_ids) 

        # Sanity check: ensure there are no duplicates in the group number associated with the primary_gal_ids
        primary_gal_ids = np.int32(primary_gal_ids)
        grnr_arr = self.cat.SubhaloGrNr[primary_gal_ids]
        if np.size(np.unique(grnr_arr)) < np.size(grnr_arr): print "had non-uniques!"

        return grnr_arr


    def read_data(self,fof_num,sub_num=-1,parttype=0,verbose=False):
        data_dict = {}

        if (sub_num>=0) & (fof_num < 0):
            off = self.HaloOffset[sub_num, parttype]
            left = self.cat.SubhaloLenType[sub_num, parttype]
            # if (verbose):
            #     print "READHALO: nr / particle # / mass :", sub_num, self.cat.SubhaloLenType[sub_num, parttype], self.cat.SubhaloMassType[sub_num, parttype].astype("float64")
        if (fof_num>=0) & (sub_num < 0):
            off = self.GroupOffset[fof_num, parttype]
            left = self.cat.GroupLenType[fof_num, parttype]
            # if (verbose):
            #     print "READHALO: nr / particle # / mass :", fof_num, self.cat.GroupLenType[fof_num, parttype], self.cat.GroupMassType[fof_num, parttype].astype("float64")
        if (sub_num>=0) & (fof_num>=0):
            real_sub_num = sub_num + self.cat.GroupFirstSub[fof_num]
            off = self.HaloOffset[real_sub_num, parttype]
            left = self.cat.SubhaloLenType[real_sub_num, parttype]
            if (verbose):
                    print "READHALO: nr / particle # / mass :", real_sub_num, self.cat.SubhaloLenType[real_sub_num, parttype], self.cat.SubhaloMassType[real_sub_num, parttype].astype("float64")
        
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


    def save_data(self,grnr,data_dict,verbose=False):
        filepath = self.savebase+"s{}/{}.hdf5".format(self.snapnum,str(int(grnr)).zfill(5))

        if verbose:
            print "Saving CGM snapshot for grnr {} to {}".format(grnr,filepath)

        f=h5py.File(filepath,'a')

        grp = f.create_group("Header")
        grp.attrs["grp_id"] = grnr
        grp.attrs["grp_mass"] = self.cat.Group_M_Crit200[grnr]
        grp.attrs["grp_pos"] = self.cat.GroupPos[grnr]
        grp.attrs["grp_Rvir"] = self.cat.Group_R_Crit200[grnr]

        p_grp = f.create_group('PartType0')
        for key in data_dict:
            p_grp.create_dataset(key,data=data_dict[key])
        f.close()



if __name__ == '__main__':
    illustris_CGM_extract()
