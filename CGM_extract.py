import numpy as np
from mpi4py import MPI
import h5py
import os.path
from sys import getsizeof


from scipy.spatial import cKDTree

import readsubfHDF5

from units import AREPO_units
AU = AREPO_units()


class CGM_extract_class:
    def __init__(self,overwrite=False):
        # self.run = "c0_512"
        # self.fnummax = 32
        # self.base="/n/hernquistfs1/Illustris/SmallBox/GFM/Production/Cosmo/Cosmo0_V6/L25n512/output/"
        # self.snapnum_arr = np.array([120])
        self.very_large_sim = True

        self.run = "ill1"
        self.fnummax=512 # future: find this automatically instead of having it as an input
        self.base="/n/ghernquist/Illustris/Runs/Illustris-1/output/"
        self.snapnum_arr = np.array([120])
        
        group_min_mass = 10.**11
        group_max_mass = 10.**18.
        self.dat_str_list = ["Masses","Coordinates","GFM_Metals","GFM_Metallicity","Velocities","Density","Volume","InternalEnergy","ElectronAbundance","NeutralHydrogenAbundance","SmoothingLength"] #future: add "Radius"
        savebase = '/n/home04/jsuresh/scratch1/Test/'
        # savebase = '/n/home04/jsuresh/data1/Projects/Feedback_and_CGM/CGM_new/data/CGM_snaps/'

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
                CGM_snapdir = savebase+"{}/s{}/".format(self.run,snapnum)
                if not os.path.isdir(CGM_snapdir):
                    print "Trying to create {}".format(CGM_snapdir)
                    os.mkdir(CGM_snapdir)

                # Get header information before proceeding:
                fname = self.base+"snapdir_"+str(snapnum).zfill(3)+"/snap_"+str(snapnum).zfill(3)+".0.hdf5"
                print "fname ",fname
                f = h5py.File(fname,'r')
                redshift=f["Header"].attrs["Redshift"]
                hubble=f["Header"].attrs["HubbleParam"]
                box=f["Header"].attrs["BoxSize"]
                #print "box: ",box
                omegam=f["Header"].attrs["Omega0"]
                omegal=f["Header"].attrs["OmegaLambda"]
                #print "redshift: ", redshift
                f.close()

                cat=readsubfHDF5.subfind_catalog(self.base,snapnum,subcat=False,keysel=["Group_M_Crit200","GroupPos","Group_R_Crit200"]) #,"GroupBHMass","GroupBHMdot"

                # Select by minimum group mass threshold:
                grp_mass = np.array(cat.Group_M_Crit200)
                m = AU.PhysicalMass(grp_mass)
                mass_select = np.logical_and(m > group_min_mass, m < group_max_mass)

                grp_mass = grp_mass[mass_select]
                grp_ids = np.arange(np.float64(cat.ngroups))
                grp_ids = grp_ids[mass_select]
                grp_pos = np.array(cat.GroupPos)[mass_select]
                grp_Rvir = np.array(cat.Group_R_Crit200)[mass_select]
                # grp_BHmass = np.array(cat.GroupBHMass)[mass_select]
                # grp_BHMdot = np.array(cat.GroupBHMdot)[mass_select]
                n_selected_groups = np.float32(np.size(grp_mass))

                if not overwrite:
                    # First remove all groups which already have data files output:
                    keep = np.ones_like(grp_ids,dtype=bool)
                    for i in np.arange(n_selected_groups):
                        grp_id = grp_ids[i]
                        filepath = savebase+ "{}/s{}/{}.hdf5".format(self.run,snapnum,str(int(grp_id)).zfill(5))
                        if os.path.isfile(filepath):
                            if self.rank == 0:
                                print "File {} already exists!  Skipping...".format(filepath)
                                keep[i] = False
                            else:
                                pass
                    grp_ids = grp_ids[keep]
                    grp_pos = grp_pos[keep]
                    grp_Rvir = grp_Rvir[keep]
                    grp_mass = grp_mass[keep]
                    n_selected_groups = np.float32(np.size(grp_mass))

                print "CGM snapshots will be written for the following group ids: {}".format(grp_ids)

            else:
                redshift = None
                hubble = None
                box = None
                omegam = None
                omegal = None
                grp_mass = None
                m = None
                grp_ids = None
                grp_pos = None
                grp_Rvir = None
                n_selected_groups = None

            if self.size > 1:
                # Broadcast necessary data from root process to other processes
                redshift = comm.bcast(redshift,root=0)
                hubble = comm.bcast(hubble,root=0)
                box = comm.bcast(box,root=0)
                omegam = comm.bcast(omegam,root=0)
                omegal = comm.bcast(omegal,root=0)
                grp_mass = comm.bcast(grp_mass,root=0)
                m = comm.bcast(m,root=0)
                grp_ids = comm.bcast(grp_ids,root=0)
                grp_pos = comm.bcast(grp_pos,root=0)
                grp_Rvir = comm.bcast(grp_Rvir,root=0)
                n_selected_groups = comm.bcast(n_selected_groups,root=0)



            if self.very_large_sim:
                # Read in positions from entire snapshot:
                pos = self.subprocess_load_data(snapnum,"Coordinates")
                if self.rank == 0: print "Done loading full position array in snapshot"

                # Build KDtree of gas positions:
                gas_kdtree = cKDTree(pos)
                print "Rank {} is done building its gas KD tree".format(self.rank)

                # Loop over all group positions, and find indices of gas close to this group position
                index_dict = {}
                for i in np.arange(n_selected_groups):
                    grp_id = grp_ids[i]
                    gpos  = grp_pos[i]
                    # r200  = grp_Rvir[i]

                    dist_thresh = AU.CodePosition(500.*np.sqrt(3),redshift)
                    ind = self.get_gas_ind(gpos,dist_thresh,gas_kdtree,box)
                    index_dict[grp_id] = np.copy(ind)

                # Now we have full index dictionary for all groups.  We no longer need the KD-tree.
                try: print "Size of index dictionary in bytes: ",sys.getsizeof(index_dict)
                except: print "try1 failed"
                try: print "Size of KD-tree in bytes: ",sys.getsizeof(gas_kdtree)
                except: print "try2 failed"
                del gas_kdtree

                # Create data structure which will house data for all groups simultaneously:
                group_dict = {}
                for i in np.arange(n_selected_groups):
                    grp_id = grp_ids[i]
                    # print "creating dict for grpid ",grp_id
                    group_dict[grp_id] = {}

                # Now get position information for every group, since we already have the position array in memory.
                # Then delete the pos array, since it is probably very large.
                for i in np.arange(n_selected_groups):
                    grp_id = grp_ids[i]
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
                        for i in np.arange(n_selected_groups):
                            grp_id = grp_ids[i]
                            ind = index_dict[grp_id]
                            group_dict[grp_id][dat_str] = dat[ind]
                        del dat
                    print "Rank {} is done with the {} field for all groups!".format(self.rank,dat_str)


                # We have now loaded all of the data for all of the groups!  Time to send it all to the rank-0 process.
                if self.rank == 0:
                    save_dict = group_dict
                if self.size > 1:
                    comm.Barrier()
                    print "Rank {} hit the comm barrier!".format(self.rank)
                    for temp_rank in range(1,self.size):
                        if self.rank == temp_rank:
                            comm.send(group_dict,dest=0,tag=self.rank)
                        elif self.rank == 0:
                            hold = comm.recv(source=temp_rank,tag=self.rank)
                            # Now combine "hold" with the dictionary we already have

                            # Future: Hide the following loop in a function: self.append_subprocess_data(save_dict,hold)
                            # For each group, check whether it is in 'hold'
                            for i in np.arange(n_selected_groups):
                                grp_id = grp_ids[i]
                                if hold.has_key(grp_id):
                                    # For each data field, check whether it is in the group data
                                    for dat_str in self.dat_str_list:
                                        if hold[grp_id].has_key(dat_str):
                                            # Then we want to add this to the dictionary ("save_dict") we already have:
                                            if save_dict.has_key(grp_id):
                                                if save_dict[grp_id].has_key(dat_str):
                                                    save_dict[grp_id] = np.append(save_dict[grp_id][dat_str],np.copy(hold[grp_id][dat_str]),axis=0)
                                                else:
                                                    save_dict[grp_id] = np.copy(hold[grp_id][dat_str])
                                            else:
                                                save_dict[grp_id] = hold[grp_id]
                                        else:
                                            pass
                                else:
                                    pass


                # Now data for this halo has been sent to the rank-0 process.  Time to save it.
                if self.rank==0:
                    print "grp_ids ",grp_ids
                    for i in np.arange(n_selected_groups):
                        grp_id = grp_ids[i]
                        filepath = savebase+ "{}/s{}/{}.hdf5".format(self.run,snapnum,str(int(grp_id)).zfill(5))
                        print "Saving group file now to: {}".format(filepath)
                        f=h5py.File(filepath,'a')

                        grp = f.create_group("Header")
                        grp.attrs["hubble"]=hubble
                        grp.attrs["omegam"]=omegam
                        grp.attrs["omegal"]=omegal
                        grp.attrs["redshift"]=redshift
                        grp.attrs["box"]=box

                        grp.attrs["grp_id"] = grp_id
                        grp.attrs["grp_mass"] = grp_mass[i]
                        grp.attrs["grp_pos"] = grp_pos[i]
                        grp.attrs["grp_Rvir"] = grp_Rvir[i]
                        # grp.attrs["grp_BHmass"] = grp_BHmass[i]
                        # grp.attrs["grp_BHMdot"] = grp_BHMdot[i]

                        p_grp = f.create_group('PartType0')
                        for dat_str in self.dat_str_list:
                            p_grp.create_dataset(dat_str,data=save_dict[grp_id][dat_str])
                        f.close()


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
                    grp_id = grp_ids[i]
                    gpos  = grp_pos[i]
                    # r200  = grp_Rvir[i]

                    # Check whether CGM data file already exists:
                    filepath = savebase+ "{}/s{}/{}.hdf5".format(run,snapnum,str(int(grp_id)).zfill(5))
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
                            grp.attrs["grp_mass"] = grp_mass[i]
                            grp.attrs["grp_pos"] = grp_pos[i]
                            grp.attrs["grp_Rvir"] = grp_Rvir[i]
                            grp.attrs["grp_BHmass"] = grp_BHmass[i]
                            grp.attrs["grp_BHMdot"] = grp_BHMdot[i]

                            p_grp = f.create_group('PartType0')
                            for key in save_dict:
                                #print save_dict[key]
                                p_grp.create_dataset(key,data=save_dict[key])
                            f.close()

    #====================================================================================================
    def get_gas_ind(self,grp_pos,dist_thresh,gas_kdtree,box):
        """ Extracts the indices of the gas which are located within dist_thresh of a group position """
        ind = gas_kdtree.query_ball_point(grp_pos, dist_thresh)

        #print "dist_thresh ",dist_thresh
        #print "ind ",ind
        # Check whether this group is close to any of the walls of the box.
        if grp_pos[0] > box-dist_thresh:
            bc_grp_pos = grp_pos - np.array([box,0,0])
            ind_bc = gas_kdtree.query_ball_point(bc_grp_pos,dist_thresh)
            ind = np.append(ind,ind_bc)

        if grp_pos[0] < dist_thresh:
            bc_grp_pos = grp_pos + np.array([box,0,0])
            ind_bc = gas_kdtree.query_ball_point(bc_grp_pos,dist_thresh)
            ind = np.append(ind,ind_bc)

        if grp_pos[1] > box-dist_thresh:
            bc_grp_pos = grp_pos - np.array([0,box,0])
            ind_bc = gas_kdtree.query_ball_point(bc_grp_pos,dist_thresh)
            ind = np.append(ind,ind_bc)

        if grp_pos[1] < dist_thresh:
            bc_grp_pos = grp_pos + np.array([0,box,0])
            ind_bc = gas_kdtree.query_ball_point(bc_grp_pos,dist_thresh)
            ind = np.append(ind,ind_bc)

        if grp_pos[2] > box-dist_thresh:
            bc_grp_pos = grp_pos - np.array([0,0,box])
            ind_bc = gas_kdtree.query_ball_point(bc_grp_pos,dist_thresh)
            ind = np.append(ind,ind_bc)

        if grp_pos[2] < dist_thresh:
            bc_grp_pos = grp_pos + np.array([0,0,box])
            ind_bc = gas_kdtree.query_ball_point(bc_grp_pos,dist_thresh)
            ind = np.append(ind,ind_bc)

        #n_found = np.size(ind)
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
    # def append_subprocess_data(self,save_dict,hold):



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
