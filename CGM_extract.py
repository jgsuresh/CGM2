import numpy as np
#import sys
from mpi4py import MPI
import h5py
import os.path

from scipy.spatial import cKDTree

import readsubfHDF5

from units import AREPO_units
AU = AREPO_units()


def main(run=None,base=None):
    # run = "ill1"
    run = "no_metal_cooling"
    base="/n/hernquistfs1/jsuresh/Runs/{}/output/".format(run)
    fnummax=8
    #snapnum_arr = np.array([54,60,68,85,114,135])
    #snapnum_arr = np.array([1,3,5,7])
    snapnum_arr = np.array([5])
    # snapnum_arr = np.array([60])
    group_min_mass = 10.**11.8 #10.**12.2
    group_max_mass = 10.**12.2
    #dat_str_list = ["Masses","Coordinates","GFM_Metals","GFM_Metallicity","Velocities","Density","Volume","InternalEnergy","ElectronAbundance","Radius","NeutralHydrogenAbundance"]
    dat_str_list = ["Masses","Coordinates","GFM_Metals","GFM_Metallicity","Velocities","Density","Volume","InternalEnergy","ElectronAbundance","NeutralHydrogenAbundance","SmoothingLength"]
    savebase = '/n/home04/jsuresh/CGM_new/data/CGM_snaps/'


    if run == "gam_25_BH":
        base="/n/hernquistfs1/jsuresh/Runs/gam_25_BH/output/"
        fnummax=8

    if run == "c0_128":
        base="/n/hernquistfs1/mvogelsberger/projects/GFM/Production/Cosmo/Cosmo0_V6/L25n128/output/"
        fnummax=8
        #boxsize = 25000
    elif run == "c0_256":
        base="/n/hernquistfs1/mvogelsberger/projects/GFM/Production/Cosmo/Cosmo0_V6/L25n256/output/"
        fnummax=8
    elif run == "c0_fw_256":
        base="/n/hernquistfs1/mvogelsberger/projects/GFM/Production/Cosmo/Cosmo0_V6_fastWinds/L25n256/output/"
        fnummax=8
    elif run == "c0_sw_256":
        base="/n/hernquistfs1/mvogelsberger/projects/GFM/Production/Cosmo/Cosmo0_V6_strongWinds/L25n256/output/"
        fnummax=8
        #boxsize = 25000
    elif run == "c2_256":
        base="/n/hernquistfs1/mvogelsberger/projects/GFM/Production/Cosmo/Cosmo2_V6/L25n256/output/"
        fnummax=8
        #boxsize = 25000
    elif run == "c0_512":
        base="/n/hernquistfs1/mvogelsberger/projects/GFM/Production/Cosmo/Cosmo0_V6/L25n512/output/"
        fnummax=32
    elif run == "c4_512":
        base="/n/hernquistfs1/spb/Cosmo/Cosmo4_V6/L25n512/output/"
        fnummax=8
    elif run == "c3_256":
        base="/n/hernquistfs1/mvogelsberger/projects/GFM/Production/Cosmo/Cosmo3_V6/L25n256/output/"
        fnummax=8
    elif run == "c3_512":
        base="/n/hernquistfs1/spb/Cosmo/Cosmo3_V6/L25n512/output/"
        fnummax=8
    elif run == "c5_256":
        base="/n/home04/jsuresh/runs/Cosmo5_V6/output/"
        fnummax=8
    if run == 'ill3':
        base="/n/ghernquist/Illustris/Runs/Illustris-3/output/"
        fnummax=32
        #boxsize = 75000
    elif run == 'ill2':
        base="/n/ghernquist/Illustris/Runs/Illustris-2/output/"
        fnummax=256
        #boxsize = 75000
    elif run == 'ill1':
        base="/n/ghernquist/Illustris/Runs/Illustris-1/output/"
        fnummax=512
        #boxsize = 75000


    # Here's where the magic happens
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    print "my rank = {}".format(rank)
    size = comm.Get_size()
    print "my size = {}".format(size)
    #print "done with MPI comm/rank/size initialization!"

    if fnummax % size != 0:
        raise Exception("# of processes does not divide into # of subfiles!")


    for snapnum in snapnum_arr:
        if rank == 0:
            # Create necessary folder if it does not exist:
            CGM_snapdir = savebase+"{}/s{}/".format(run,snapnum)
            if not os.path.isdir(CGM_snapdir):
                print "Trying to create {}".format(CGM_snapdir)
                os.mkdir(CGM_snapdir)

            # Get header information before proceeding:
            fname = base+"snapdir_"+str(snapnum).zfill(3)+"/snap_"+str(snapnum).zfill(3)+".0.hdf5"
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

            cat=readsubfHDF5.subfind_catalog(base,snapnum,subcat=False)
            #cat=readsubfHDF5.subfind_catalog(base,135,keysel=["GroupMass","GroupPos","Group_R_Crit200"])

            # Select by minimum group mass threshold:
            grp_mass = np.array(cat.Group_M_Crit200) #check what hubble is
            m = AU.PhysicalMass(grp_mass)
            #mass_select = (m > group_min_mass)
            mass_select = np.logical_and(m > group_min_mass, m < group_max_mass)

            grp_mass = grp_mass[mass_select]
            grp_ids = np.arange(np.float64(cat.ngroups))
            grp_ids = grp_ids[mass_select]
            grp_pos = np.array(cat.GroupPos)[mass_select]
            grp_Rvir = np.array(cat.Group_R_Crit200)[mass_select]
            #grp_BHmass = np.array(cat.GroupBHMass)[mass_select]
            #grp_BHMdot = np.array(cat.GroupBHMdot)[mass_select]
            n_selected_groups = np.float32(np.size(grp_mass))
            del cat
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

        if size > 1:
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




        # Read in all data from snapshot a single time:
        local_dict = {}
        for fnum in range(0,fnummax):
            if (fnum % size == rank):
                print "Projecting file / task:", fnum, rank
                fname = base+"snapdir_"+str(snapnum).zfill(3)+"/snap_"+str(snapnum).zfill(3)+"."+str(fnum)+".hdf5"

                print "Reading from {}".format(fname)
                f = h5py.File(fname,'r')
                for dat_str in dat_str_list:
                    dat = np.array(f['PartType0'][dat_str],dtype=np.float32)

                    if local_dict.has_key(dat_str): 
                        local_dict[dat_str] = np.append(local_dict[dat_str],dat,axis=0)
                    else:
                        local_dict[dat_str] = np.array(f['PartType0'][dat_str],dtype=np.float32)

        if rank == 0:
            print "Done loading full snapshot"

        # Build KDtree of gas positions:
        pos = local_dict['Coordinates']
        gas_kdtree = cKDTree(pos)

        # Loop over all group positions, and find indices of gas close to this group position
        for i in np.arange(n_selected_groups):
            grp_id = grp_ids[i]
            gpos  = grp_pos[i]
            r200  = grp_Rvir[i]

            # Check whether file exists:
            filepath = savebase+ "{}/s{}/{}.hdf5".format(run,snapnum,str(int(grp_id)).zfill(5))
            if os.path.isfile(filepath):
                if rank == 0:
                    print "File {} already exists!  Skipping...".format(filepath)
                else:
                    pass
            else:
                dist_thresh = AU.CodePosition(500.*np.sqrt(3),redshift) #3*r200
                ind = get_gas_ind(gpos,dist_thresh,gas_kdtree,box)

                send_dict = {}
                if np.size(ind) > 0:
                    for dat_str in dat_str_list:
                        send_dict[dat_str] = local_dict[dat_str][ind]

                # Now that data has been compiled for this process, send it to rank-0 process
                if rank == 0:
                    save_dict = send_dict
                if size > 1:
                    comm.Barrier()
                    print "comm barrier!"
                    for dat_str in dat_str_list:
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
                    #grp.attrs["grp_BHmass"] = grp_BHmass[i]
                    #grp.attrs["grp_BHMdot"] = grp_BHMdot[i]

                    p_grp = f.create_group('PartType0')
                    for key in save_dict:
                        #print save_dict[key]
                        p_grp.create_dataset(key,data=save_dict[key])
                    f.close()




def get_gas_ind(grp_pos,dist_thresh,gas_kdtree,box):
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


if __name__ == '__main__':
    main()
