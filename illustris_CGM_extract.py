# Unlike CGM_extract, which extracts all gas within a specific radius of each desired halo, this
# script focuses on extracting gas that is bound to the FoF groups, thus exploiting the natural
# format of the snapshots themselves.

import numpy as np
#import sys
from mpi4py import MPI
import h5py
import os.path

from scipy.spatial import cKDTree

import readsubfHDF5

from units import AREPO_units
AU = AREPO_units()

import readhaloHDF5


class illustris_CGM_extract:
    def __init__(self):
        run = "ill1"
        self.snapnum = 120

        # group_min_mass = 10.**11.8 
        # group_max_mass = 10.**12.2
        # gal_min_stellar_mass = 10.**10.
        # gal_max_stellar_mass = 10.**11.6
        dat_str_list = ["Masses","Coordinates","GFM_Metals","GFM_Metallicity","Velocities","Density","Volume","InternalEnergy","ElectronAbundance","NeutralHydrogenAbundance","SmoothingLength"]
        self.block_list = ["MASS","POS ","GMET","GZ  ","VEL ","RHO ","VOL ","U   ","NE  ","HSML"]
        self.savebase = '/n/home04/jsuresh/scratch1/AREPOfest/data/CGM_snaps/'

        if run == 'ill3': self.base="/n/ghernquist/Illustris/Runs/Illustris-3/output/"
        elif run == 'ill2': self.base="/n/ghernquist/Illustris/Runs/Illustris-2/output/"
        elif run == 'ill1': self.base="/n/ghernquist/Illustris/Runs/Illustris-1/output/"

        self.setup_tables()
        self.read_data(13,verbose=True)

        # Let's try using readhaloHDF5.  I'm curious about how fast it is!
        # dat_str_list = ["Masses","Coordinates","GFM_Metals","GFM_Metallicity","Velocities","Density","Volume","InternalEnergy","ElectronAbundance","NeutralHydrogenAbundance","SmoothingLength"]
        # type=4
        # grpnr=13
        # subnr=-1 #extract full FoF group
        # readhaloHDF5.reset()
        # pos=readhaloHDF5.readhalo(base, "snap", snapnum, "POS ", type, grpnr, subnr, long_ids=True, double_output=False)



    def setup_tables(self,parttype=0,verbose=False,double_output=False):
        # global FlagRead, cat, GroupOffset, HaloOffset, multiple, filename, Parttype, FileTypeNumbers, FileNum
 
        if (verbose):
            print "READHALO: INITIAL READ"
            print "READHALO: Parttype = ", parttype

        #read in catalog
        self.cat = readsubfHDF5.subfind_catalog(base, self.snapnum, long_ids=long_ids, double_output=double_output, keysel=["GroupLenType","GroupNsubs","GroupFirstSub","SubhaloLenType","SubhaloMassType"])

        if (self.cat.ngroups==0):
            if (verbose):
                print "READHALO: no groups in catalog... returning"
            return

        
        if (FlagRead==False):
            self.GroupOffset = np.zeros([self.cat.ngroups, 6], dtype="int64")
            self.HaloOffset = np.zeros([self.cat.nsubs, 6], dtype="int64")
            filename = base+"/"+self.snapbase+"_"+str(self.snapnum).zfill(3)
            self.multiple=False
            if (os.path.exists(filename+".hdf5")==False):
                    filename = base+"/snapdir_"+str(self.snapnum).zfill(3)+"/"+self.snapbase+"_"+str(self.snapnum).zfill(3)+"."+str(0)
                    self.multiple=True
            if (os.path.exists(filename+".hdf5")==False):
                    print "READHALO: [error] file not found : ", filename
                    sys.exit()

            FlagRead=True


        #construct offset tables
        k=0
        for i in range(0, self.cat.ngroups):
            if (i>0):
                self.GroupOffset[i, parttype] =  self.GroupOffset[i-1, parttype] + cat.GroupLenType[i-1, parttype]
            if (self.cat.GroupNsubs[i]>0):
                self.HaloOffset[k, parttype] = self.GroupOffset[i, parttype]
                k+=1
                for j in range(1, cat.GroupNsubs[i]):
                    self.HaloOffset[k, parttype] =  self.HaloOffset[k-1, parttype] + cat.SubhaloLenType[k-1, parttype]
                    k+=1
        if (k!=self.cat.nsubs):
            print "READHALO: problem with offset table", k, cat.nsubs
            sys.exit()

        #construct file tables
        if (self.multiple):
            filename = base+"/snapdir_"+str(self.snapnum).zfill(3)+"/"+self.snapbase+"_"+str(self.snapnum).zfill(3)+"."+str(0)
        else:
            filename = base+"/"+self.snapbase+"_"+str(self.snapnum).zfill(3)

        head = snapHDF5.snapshot_header(filename)
        self.FileNum = head.filenum

        self.FileTypeNumbers = np.zeros([self.FileNum, 6], dtype="int64") 
        cumcount = np.zeros(6, dtype="int64")

        for fnr in range(0, self.FileNum-1):
            if (self.multiple):
                filename = base+"/snapdir_"+str(self.snapnum).zfill(3)+"/"+self.snapbase+"_"+str(self.snapnum).zfill(3)+"."+str(fnr)
            else:
                filename = base+"/"+self.snapbase+"_"+str(self.snapnum).zfill(3)

            if (verbose):
                print "READHALO: initial reading file :", filename

            head = snapHDF5.snapshot_header(filename)
    
            cumcount[:] += head.npart[:]
            self.FileTypeNumbers[fnr+1, :] = cumcount[:]


    def read_data(self,fof_num,sub_num=-1,parttype=0):
        data_dict = {}

        if (sub_num>=0) & (fof_num < 0):
            off = self.HaloOffset[sub_num, parttype]
            left = self.cat.SubhaloLenType[sub_num, parttype]
            if (verbose):
                print "READHALO: nr / particle # / mass :", sub_num, self.cat.SubhaloLenType[sub_num, parttype], self.cat.SubhaloMassType[sub_num, parttype].astype("float64")
        if (fof_num>=0) & (sub_num < 0):
            off = self.GroupOffset[fof_num, parttype]
            left = self.cat.GroupLenType[fof_num, parttype]
            if (verbose):
                print "READHALO: nr / particle # / mass :", fof_num, self.cat.GroupLenType[fof_num, parttype], self.cat.GroupMassType[fof_num, parttype].astype("float64")
        if (sub_num>=0) & (fof_num>=0):
            real_sub_num = sub_num + self.cat.GroupFirstSub[fof_num]
            off = self.HaloOffset[real_sub_num, parttype]
            left = cat.SubhaloLenType[real_sub_num, parttype]
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
                filename = base+"/snapdir_"+str(self.snapnum).zfill(3)+"/"+self.snapbase+"_"+str(self.snapnum).zfill(3)+"."+str(fnr)
            else:
                filename = base+"/"+self.snapbase+"_"+str(self.snapnum).zfill(3)
        
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
                        if (verbose):
                            print "Reading block again: ",block_name
                        data_dict[block_name] = np.append(data_dict[block_name], snapHDF5.read_block(filename, block_name, parttype, slab_start=start, slab_len=count), axis=0)
                    left -= count
                    off += count
            if (left==0):
                break
            off -= nloc

        return data




if __name__ == '__main__':
    main()
