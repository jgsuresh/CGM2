import os
import time
import numpy as np
import h5py

import matplotlib
matplotlib.use('PDF')
from matplotlib import rc
#rc('text', usetex=True)
rc('font', family='serif')
import matplotlib.pyplot as plt
import matplotlib.ticker
import brewer2mpl
from itertools import cycle

from units import AREPO_units
AU = AREPO_units()
import readsubfHDF5

# import SphMap as sm 

import glob

snapnum = 68

fig_base = '/n/home04/jsuresh/CGM_new/data/figs/'
npz_base = '/n/home04/jsuresh/CGM_new/data/npz/'

mv_base = "/n/hernquistfs1/Illustris/SmallBox/GFM/Production/Cosmo/"
# snapbase_list = [mv_base+"Cosmo0_V6/L25n256/output/",mv_base+"TurnOnOff/Cosmo0_V6_MetSteWiStrQuRad/L25n256/output/" ]
# run_list = ["c0","norad"]
snapbase_list = [mv_base+"Cosmo0_V6/L25n256/output/"]
run_list = ["c0"]

def _fix_pos(grp_pos,pos,boxsize):

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

plt.figure()

if False:
    for i in np.arange(len(snapbase_list)):
        print "snapbase_list[i] ",snapbase_list[i]
        cat = readsubfHDF5.subfind_catalog(snapbase_list[i],snapnum,keysel=['GroupPos','Group_M_Crit200','Group_R_Crit200'],subcat=False)

        grp_m = AU.PhysicalMass(cat.Group_M_Crit200)
        mass_cut = np.logical_and(grp_m >= 10.**11.8,grp_m < 10.**12.2)

        grp_ids = np.arange(np.size(grp_m))

        grp_m = grp_m[mass_cut]
        grp_Rvir = cat.Group_R_Crit200[mass_cut]
        grp_pos = cat.GroupPos[mass_cut]
        grp_ids = grp_ids[mass_cut]
        n_grps = np.size(grp_m)

        #Open snapshot and get gas positions and masses:
        # gas_pos = np.zeros([0,3])
        # gas_mass = np.zeros(0)
        for fn in np.arange(8):
            f = h5py.File(snapbase_list[i]+"snapdir_068/snap_068.{}.hdf5".format(str(fn)),'r')
            if fn == 0:
                gas_mass = np.copy(f['PartType0']['Masses'])
                gas_pos = np.copy(f['PartType0']['Coordinates'])
                star_mass = np.copy(f['PartType4']['Masses'])
                star_pos = np.copy(f['PartType4']['Coordinates'])
                star_age = np.copy(f['PartType4']['GFM_StellarFormationTime'])
            else:
                gas_mass = np.append(gas_mass,np.copy(f['PartType0']['Masses']))
                gas_pos = np.append(gas_pos,np.copy(f['PartType0']['Coordinates']),axis=0)
                star_mass = np.append(star_mass,np.copy(f['PartType4']['Masses']))
                star_pos = np.append(star_pos,np.copy(f['PartType4']['Coordinates']),axis=0)
                star_age = np.append(star_age,np.copy(f['PartType4']['GFM_StellarFormationTime']))

        grp_gas_m = np.zeros(n_grps)
        grp_star_m = np.zeros(n_grps)
        grp_wind_m = np.zeros(n_grps)
        # n_grps = 1
        for j in np.arange(n_grps):
            print "{} of {}".format(j,n_grps)
            cent = grp_pos[j]
            pos_cent = _fix_pos(cent,gas_pos,25000)
            r = AU._dist(pos_cent,np.array([0.,0.,0.]))
            r_cut = r < grp_Rvir[j]
            grp_gas_m[j] = np.sum(gas_mass[r_cut])

            pos_cent = _fix_pos(cent,star_pos,25000)
            r = AU._dist(pos_cent,np.array([0.,0.,0.]))
            r_cut = r < grp_Rvir[j]
            wind = star_age[r_cut] < 0.
            not_wind = np.logical_not(wind)
            grp_wind_m[j] = np.sum(star_mass[r_cut][wind])
            grp_star_m[j] = np.sum(star_mass[r_cut][not_wind])

        grp_gas_m = AU.PhysicalMass(grp_gas_m)
        grp_wind_m = AU.PhysicalMass(grp_wind_m)
        grp_star_m = AU.PhysicalMass(grp_star_m)
        npz_fname = "bary_abund_z2_rvir_wwind_{}.npz".format(run_list[i])

        np.savez(npz_base+npz_fname,grp_m=grp_m,grp_gas_m=grp_gas_m,grp_wind_m=grp_wind_m,grp_star_m=grp_star_m)
        # np.savez(grp_gas_m,grp_star_m,grp_mass,grp_Rvir,grp_ids)
        # scatter plot

        plt.scatter(np.log10(grp_m),grp_gas_m/grp_m,label='gas')
        plt.scatter(np.log10(grp_m),grp_star_m/grp_m,color='red',label='stars')
        plt.scatter(np.log10(grp_m),grp_wind_m/grp_m,color='green',label='wind')
        plt.scatter(np.log10(grp_m),(grp_gas_m+grp_star_m)/grp_m,color='black',label='total')
        plt.savefig(fig_base+"wind_baryon_check.pdf")


if True:
    # f_c0 = np.load(npz_base+'bary_abund_z2_3rvir_c0.npz')
    # grp_m = f_c0['grp_m']
    # grp_gas_m = f_c0['grp_gas_m']
    # grp_star_m = f_c0['grp_star_m']
    # plt.scatter(np.log10(grp_m),(grp_gas_m+grp_star_m)/grp_m,color='blue',label='c0')

    # f_c0 = np.load(npz_base+'bary_abund_z2_3rvir_norad.npz')
    # grp_m = f_c0['grp_m']
    # grp_gas_m = f_c0['grp_gas_m']
    # grp_star_m = f_c0['grp_star_m']
    # plt.scatter(np.log10(grp_m),(grp_gas_m+grp_star_m)/grp_m,color='red',label='no radio')
    # plt.savefig(fig_base+"halogas_z2_radio_comp_3rvir.pdf")

    f_c0 = np.load(npz_base+'bary_abund_z2_rvir_wwind_c0.npz')
    grp_m = f_c0['grp_m']
    grp_gas_m = f_c0['grp_gas_m']
    grp_star_m = f_c0['grp_star_m']
    grp_wind_m = f_c0['grp_wind_m']
    plt.scatter(np.log10(grp_m),np.log10(grp_gas_m),color='blue',label='gas')
    # plt.scatter(np.log10(grp_m),np.log10(grp_star_m),color='red',label='stars')
    plt.scatter(np.log10(grp_m),np.log10(grp_wind_m),color='green',label='wind')
    plt.scatter(np.log10(grp_m),np.log10((grp_gas_m+grp_star_m+grp_wind_m)),color='black',label='total')
    plt.legend()
    plt.savefig(fig_base+"wind_baryon_check.pdf")












