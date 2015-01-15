import numpy as np
from mpi4py import MPI
# import mpi4py.rc
# mpi4py.rc.threaded = False 
import h5py
import os.path
from sys import getsizeof


from scipy.spatial import cKDTree

import readsubfHDF5

from units import AREPO_units
AU = AREPO_units()


class CGM_extract_class:
    def __init__(self,overwrite=False):
        # self.run = "c0_128"
        # self.fnummax = 8
        # self.base="/n/hernquistfs1/Illustris/SmallBox/GFM/Production/Cosmo/Cosmo0_V6/L25n128/output/"
        # self.snapnum_arr = np.array([120])
        self.very_large_sim = True

        self.run = "ill1"
        self.fnummax=512 # future: find this automatically instead of having it as an input
        self.base="/n/ghernquist/Illustris/Runs/Illustris-1/output/"
        self.snapnum_arr = np.array([120])
        
        group_min_mass = 10.**11
        group_max_mass = 10.**18.
        self.dat_str_list = ["Masses","Coordinates","GFM_Metals","GFM_Metallicity","Velocities","Density","Volume","InternalEnergy","ElectronAbundance","NeutralHydrogenAbundance","SmoothingLength"] #future: add "Radius"
        self.savebase = '/n/home04/jsuresh/scratch1/Test/'
        # self.savebase = '/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/CGM_new/data/CGM_snaps/'

        # Here's where the magic happens
        comm = MPI.COMM_WORLD
        self.rank = comm.Get_rank()
        print "my rank = {}".format(self.rank)
        self.size = comm.Get_size()
        # print "my size = {}".format(size)
        if self.rank == 0: print "Done with MPI comm/rank/size initialization!"

        if self.fnummax % self.size != 0:
            raise Exception("# of processes does not divide into # of subfiles!")

        for snapnum in self.snapnum_arr:
            if self.rank == 0:
                # Create necessary folder if it does not exist:
                CGM_snapdir = self.savebase+"{}/s{}/".format(self.run,snapnum)
                if not os.path.isdir(CGM_snapdir):
                    print "Trying to create {}".format(CGM_snapdir)
                    os.mkdir(CGM_snapdir)

                # Get header information before proceeding:
                fname = self.base+"snapdir_"+str(snapnum).zfill(3)+"/snap_"+str(snapnum).zfill(3)+".0.hdf5"
                print "fname ",fname
                f = h5py.File(fname,'r')
                self.redshift=f["Header"].attrs["Redshift"]
                self.hubble=f["Header"].attrs["HubbleParam"]
                self.box=f["Header"].attrs["BoxSize"]
                self.omegam=f["Header"].attrs["Omega0"]
                self.omegal=f["Header"].attrs["OmegaLambda"]
                f.close()

                cat=readsubfHDF5.subfind_catalog(self.base,snapnum,subcat=False,keysel=["Group_M_Crit200","GroupPos","Group_R_Crit200"]) #,"GroupBHMass","GroupBHMdot"

                # Select by minimum group mass threshold:
                self.grp_mass = np.array(cat.Group_M_Crit200)
                m = AU.PhysicalMass(self.grp_mass)
                mass_select = np.logical_and(m > group_min_mass, m < group_max_mass)

                self.grp_mass = self.grp_mass[mass_select]
                self.grp_ids = np.arange(np.float64(cat.ngroups))
                self.grp_ids = self.grp_ids[mass_select]
                self.grp_pos = np.array(cat.GroupPos)[mass_select]
                self.grp_Rvir = np.array(cat.Group_R_Crit200)[mass_select]
                # self.grp_BHmass = np.array(cat.GroupBHMass)[mass_select]
                # self.grp_BHMdot = np.array(cat.GroupBHMdot)[mass_select]
                self.n_selected_groups = np.float32(np.size(self.grp_mass))

                if not overwrite:
                    # First remove all groups which already have data files output:
                    keep = np.ones_like(self.grp_ids,dtype=bool)
                    for i in np.arange(self.n_selected_groups):
                        grp_id = self.grp_ids[i]
                        filepath = self.savebase+ "{}/s{}/{}.hdf5".format(self.run,snapnum,str(int(grp_id)).zfill(5))
                        if os.path.isfile(filepath):
                            if self.rank == 0:
                                print "File {} already exists!  Skipping...".format(filepath)
                                keep[i] = False
                            else:
                                pass
                    self.grp_ids = self.grp_ids[keep]
                    self.grp_pos = self.grp_pos[keep]
                    self.grp_Rvir = self.grp_Rvir[keep]
                    self.grp_mass = self.grp_mass[keep]
                    self.n_selected_groups = np.float32(np.size(self.grp_mass))

                print "CGM snapshots will be written for the following group ids: {}".format(self.grp_ids)

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

            if self.size > 1:
                # Broadcast necessary data from root process to other processes
                self.redshift = comm.bcast(self.redshift,root=0)
                self.hubble = comm.bcast(self.hubble,root=0)
                self.box = comm.bcast(self.box,root=0)
                self.omegam = comm.bcast(self.omegam,root=0)
                self.omegal = comm.bcast(self.omegal,root=0)
                self.grp_mass = comm.bcast(self.grp_mass,root=0)
                self.grp_ids = comm.bcast(self.grp_ids,root=0)
                self.grp_pos = comm.bcast(self.grp_pos,root=0)
                self.grp_Rvir = comm.bcast(self.grp_Rvir,root=0)
                self.n_selected_groups = comm.bcast(self.n_selected_groups,root=0)



            if self.very_large_sim:
                # Read in positions from entire snapshot:
                pos = self.subprocess_load_data(snapnum,"Coordinates")
                if self.rank == 0: print "Done loading full position array in snapshot"

                # Build KDtree of gas positions:
                gas_kdtree = cKDTree(pos)
                print "Rank {} is done building its gas KD tree".format(self.rank)

                # Loop over all group positions, and find indices of gas close to this group position
                index_dict = {}
                for i in np.arange(self.n_selected_groups):
                    grp_id = self.grp_ids[i]
                    gpos  = self.grp_pos[i]
                    # r200  = self.grp_Rvir[i]

                    dist_thresh = AU.CodePosition(500.*np.sqrt(3),self.redshift)
                    ind = self.get_gas_ind(gpos,dist_thresh,gas_kdtree)
                    index_dict[grp_id] = np.copy(ind)

                # Now we have full index dictionary for all groups.  We no longer need the KD-tree.
                # try: print "Size of index dictionary in bytes: ",sys.getsizeof(index_dict)
                # except: print "try1 failed"
                # try: print "Size of KD-tree in bytes: ",sys.getsizeof(gas_kdtree)
                # except: print "try2 failed"
                del gas_kdtree

                # Create data structure which will house data for all groups simultaneously:
                group_dict = {}
                for i in np.arange(self.n_selected_groups):
                    grp_id = self.grp_ids[i]
                    group_dict[grp_id] = {}

                # Now get position information for every group, since we already have the position array in memory.
                # Then delete the pos array, since it is probably very large.
                for i in np.arange(self.n_selected_groups):
                    grp_id = self.grp_ids[i]
                    ind = index_dict[grp_id]
                    group_dict[grp_id]["Coordinates"] = pos[ind]
                del pos
                print "Rank {} is done with the Coordinates field for all groups!".format(self.rank)

                # Now do the same for all other desired fields, deleting each one after we are done with it.
                for dat_str in self.dat_str_list:
                    if dat_str == "Coordinates":
                        pass
                    else:
                        dat = self.subprocess_load_data(snapnum,dat_str)
                        for i in np.arange(self.n_selected_groups):
                            grp_id = self.grp_ids[i]
                            ind = index_dict[grp_id]
                            group_dict[grp_id][dat_str] = dat[ind]
                        del dat
                    print "Rank {} is done with the {} field for all groups!".format(self.rank,dat_str)


                # We have now loaded all of the data for all of the groups!  First delete the large index_dict:
                del index_dict

                if self.size == 1:
                    # If there is only one process, then just trivially loop through all groups and save their data
                    for i in np.arange(self.n_selected_groups):
                        grp_id = self.grp_ids[i]
                        save_dict = group_dict[grp_id]
                        self.save_group_data(snapnum,i,save_dict)

                elif self.size > 1:
                    # If there are subprocesses, then send the data to the rank-0 process, where it will be saved.
                    comm.Barrier()
                    print "Rank {} hit the comm barrier!".format(self.rank)
                    # We now send data from the subprocesses to the rank-0 process, group by group.
                    for i in np.arange(self.n_selected_groups):
                        grp_id = self.grp_ids[i]
                        if self.rank == 0: 
                            save_dict = {}
                            # save_dict = group_dict[grp_id]

                        for dat_str in self.dat_str_list:
                            if group_dict.has_key(grp_id) and group_dict[grp_id].has_key(dat_str):
                                foo = np.copy(group_dict[grp_id][dat_str])
                            else:
                                foo = None
                            
                            # Now gather from all of the subprocesses to rank 0:
                            foo = comm.gather(foo,root=0)

                            if self.rank == 0:
                                savedat = self._cat_gather_list(foo)
                                save_dict[dat_str] = savedat

                        # Rank 0 now has all of the data for this group saved in save_dict.  Save to file:
                        if self.rank == 0:
                            self.save_group_data(snapnum,i,save_dict)
                            # del save_dict # HERE
                        # Now all processes remove this grp_id entry from their memory
                        # if group_dict.has_key(grp_id): 
                        #     del group_dict[grp_id]
                        # else: 
                        #     pass



            elif not self.very_large_sim:
                # Read in all data from snapshot a single time:
                local_dict = {}
                for dat_str in self.dat_str_list:
                    local_dict[dat_str] = subprocess_load_data(dat_str)
                if rank == 0: print "Done loading full snapshot"

                # Build KDtree of gas positions:
                pos = local_dict['Coordinates']
                gas_kdtree = cKDTree(pos)

                # Loop over all group positions, and find indices of gas close to this group position
                for i in np.arange(n_selected_groups):
                    grp_id = self.grp_ids[i]
                    gpos  = self.grp_pos[i]
                    # r200  = self.grp_Rvir[i]

                    # Check whether CGM data file already exists:
                    filepath = self.savebase+ "{}/s{}/{}.hdf5".format(run,snapnum,str(int(grp_id)).zfill(5))
                    if os.path.isfile(filepath):
                        if rank == 0:
                            print "File {} already exists!  Skipping...".format(filepath)
                        else:
                            pass
                    else: # If CGM data file does not exist yet, then continue:
                        dist_thresh = AU.CodePosition(500.*np.sqrt(3),redshift)
                        ind = get_gas_ind(gpos,dist_thresh,gas_kdtree,box)

                        send_dict = {}
                        if np.size(ind) > 0:
                            for dat_str in self.dat_str_list:
                                send_dict[dat_str] = local_dict[dat_str][ind]

                        # Future: hide this in another function
                        # Now that data has been compiled for this process, send it to rank-0 process
                        if rank == 0:
                            save_dict = send_dict
                        if size > 1:
                            comm.Barrier()
                            print "comm barrier!"
                            for dat_str in self.dat_str_list:
                                for temp_rank in range(1,size):
                                    if rank == temp_rank:
                                        if send_dict.has_key(dat_str):
                                            comm.send(send_dict[dat_str],dest=0,tag=i)
                                        else:
                                            comm.send('nothing',dest=0,tag=i)
                                    elif rank == 0:
                                        hold = comm.recv(source=temp_rank,tag=i)
                                        #print "this is what was received from rank {}: {}".format(temp_rank,hold)
                                        if hold != 'nothing':
                                            if save_dict.has_key(dat_str): 
                                                full_dict[dat_str] = np.append(save_dict[dat_str],hold,axis=0)
                                            else:
                                                full_dict[dat_str] = hold


                        # Now data for this halo has been sent to rank-0 process.  Time to save it.
                        if rank==0:
                            print "Saving group file now to: {}".format(filepath)
                            f=h5py.File(filepath,'a')

                            grp = f.create_group("Header")
                            grp.attrs["hubble"]=hubble
                            grp.attrs["omegam"]=omegam
                            grp.attrs["omegal"]=omegal
                            grp.attrs["redshift"]=redshift
                            grp.attrs["box"]=box

                            grp.attrs["grp_id"] = grp_id
                            grp.attrs["grp_mass"] = self.grp_mass[i]
                            grp.attrs["grp_pos"] = self.grp_pos[i]
                            grp.attrs["grp_Rvir"] = self.grp_Rvir[i]
                            grp.attrs["grp_BHmass"] = self.grp_BHmass[i]
                            grp.attrs["grp_BHMdot"] = self.grp_BHMdot[i]

                            p_grp = f.create_group('PartType0')
                            for key in save_dict:
                                #print save_dict[key]
                                p_grp.create_dataset(key,data=save_dict[key])
                            f.close()

    #====================================================================================================
    def get_gas_ind(self,gpos,dist_thresh,gas_kdtree):
        """ Extracts the indices of the gas which are located within dist_thresh of a group position """
        ind = gas_kdtree.query_ball_point(gpos, dist_thresh)

        # Check whether this group is close to any of the walls of the box.
        if gpos[0] > self.box-dist_thresh:
            bc_grp_pos = gpos - np.array([self.box,0,0])
            ind_bc = gas_kdtree.query_ball_point(bc_grp_pos,dist_thresh)
            ind = np.append(ind,ind_bc)

        if gpos[0] < dist_thresh:
            bc_grp_pos = gpos + np.array([self.box,0,0])
            ind_bc = gas_kdtree.query_ball_point(bc_grp_pos,dist_thresh)
            ind = np.append(ind,ind_bc)

        if gpos[1] > self.box-dist_thresh:
            bc_grp_pos = gpos - np.array([0,self.box,0])
            ind_bc = gas_kdtree.query_ball_point(bc_grp_pos,dist_thresh)
            ind = np.append(ind,ind_bc)

        if gpos[1] < dist_thresh:
            bc_grp_pos = gpos + np.array([0,self.box,0])
            ind_bc = gas_kdtree.query_ball_point(bc_grp_pos,dist_thresh)
            ind = np.append(ind,ind_bc)

        if gpos[2] > self.box-dist_thresh:
            bc_grp_pos = gpos - np.array([0,0,self.box])
            ind_bc = gas_kdtree.query_ball_point(bc_grp_pos,dist_thresh)
            ind = np.append(ind,ind_bc)

        if gpos[2] < dist_thresh:
            bc_grp_pos = gpos + np.array([0,0,self.box])
            ind_bc = gas_kdtree.query_ball_point(bc_grp_pos,dist_thresh)
            ind = np.append(ind,ind_bc)

        return np.uint64(ind)

    #====================================================================================================
    def subprocess_load_data(self,snapnum,dat_str):
        flag = 0
        for fnum in range(0,self.fnummax):
            if (fnum % self.size == self.rank):
                print "Task {} reading subfile {} to load {}".format(self.rank,fnum,dat_str)
                fname = self.base+"snapdir_"+str(snapnum).zfill(3)+"/snap_"+str(snapnum).zfill(3)+"."+str(fnum)+".hdf5"

                f = h5py.File(fname,'r')
                dat = np.array(f['PartType0'][dat_str],dtype=np.float32)

                if flag == 0:
                    save = dat
                    flag = 1
                elif flag == 1:
                    save = np.append(save,dat,axis=0)
        return save

    #====================================================================================================
    def _cat_gather_list(self,glist):
        # Takes list that results from comm.gather "glist" and appends the arrays inside together
        print "starting glist ",glist
        save = glist.pop(0)
        while len(glist) > 0:
            save = np.append(save,glist.pop(0),axis=0)
        return np.copy(save)


    #====================================================================================================
    def save_group_data(self,snapnum,i,save_dict):
        # Save data for this group.
        grp_id = self.grp_ids[i]

        filepath = self.savebase+ "{}/s{}/{}.hdf5".format(self.run,snapnum,str(int(grp_id)).zfill(5))
        print "Saving group file now to: {}".format(filepath)
        f=h5py.File(filepath,'a')

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
        # grp.attrs["grp_BHmass"] = self.grp_BHmass[i]
        # grp.attrs["grp_BHMdot"] = self.grp_BHMdot[i]

        p_grp = f.create_group('PartType0')
        for dat_str in self.dat_str_list:
            foo = np.array(save_dict[dat_str])
            p_grp.create_dataset(dat_str,data=foo)
        f.close()


    #====================================================================================================
    def get_gas_ind_BOX(self,grp_pos,dist_thresh,gas_pos): 
        # simple test function.  NOTE: does not handle box boundaries correctly.
        xc = np.logical_and(gas_pos[:,0] >= grp_pos[0]-dist_thresh,gas_pos[:,0] <= grp_pos[0]+dist_thresh)
        yc = np.logical_and(gas_pos[:,1] >= grp_pos[1]-dist_thresh,gas_pos[:,1] <= grp_pos[1]+dist_thresh)
        zc = np.logical_and(gas_pos[:,2] >= grp_pos[2]-dist_thresh,gas_pos[:,2] <= grp_pos[2]+dist_thresh)
        ind = np.logical_and(np.logical_and(xc,yc),zc)

        print "grp_pos[0]-dist_thresh ",grp_pos[0]-dist_thresh
        print "grp_pos[0]+dist_thresh ",grp_pos[0]+dist_thresh
        print "grp_pos[1]-dist_thresh ",grp_pos[1]-dist_thresh
        print "grp_pos[1]+dist_thresh ",grp_pos[1]+dist_thresh
        print "grp_pos[2]-dist_thresh ",grp_pos[2]-dist_thresh
        print "grp_pos[2]+dist_thresh ",grp_pos[2]+dist_thresh

        return ind


if __name__ == '__main__':
    CGM_extract_class()
