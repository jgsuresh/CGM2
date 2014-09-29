import readsubfHDF5
import h5py
import numpy as np

import matplotlib
matplotlib.use('PDF')
from matplotlib import rc
#rc('text', usetex=True)
rc('font', family='serif')
import matplotlib.pyplot as plt
import brewer2mpl

from units import AREPO_units
AU = AREPO_units()

base = "/n/ghernquist/Illustris/Runs/L75n1820FP/"
fig_base = "/n/home04/jsuresh/MgFe/data/figs/"
npz_base = "/n/home04/jsuresh/MgFe/data/npz/"

z_arr = np.array([0.,0.1,0.2,0.5,1,2,3,4])
snapnum_arr = np.array([135,127,120,103,85,68,60,54])
n_snaps = np.size(snapnum_arr)

filepath = npz_base+"MgFe_subfind_data.hdf5"
if False:
    f=h5py.File(filepath,'a')

    for snapnum in snapnum_arr:
    	print "Saving data for snapnum {}".format(snapnum)
        cat = readsubfHDF5.subfind_catalog(base,snapnum,subcat=False,keysel=['GroupMass','GroupMassType','GroupGasMetalFractions'])

        grp_mass = AU.PhysicalMass(cat.GroupMass)
        mass_cut = grp_mass > 10.**11
        grp_mass = grp_mass[mass_cut]
        gas_mass = AU.PhysicalMass(cat.GroupMassType[mass_cut,0])
        # star_mass = AU.PhysicalMass(cat.GroupMassType[mass_cut,4])
        gas_met = cat.GroupGasMetalFractions[mass_cut,:]

        #["H","He","C","N","O","Ne","Mg","Si","Fe"]
        Mg_mass = gas_mass*gas_met[:,6]
        Fe_mass = gas_mass*gas_met[:,8]

        MgFe_ratio = Mg_mass/Fe_mass

        grp = f.create_group('snap{}'.format(str(int(snapnum)).zfill(3)))
        grp.create_dataset('grp_mass',data=grp_mass)
        grp.create_dataset('gas_mass',data=gas_mass)
        grp.create_dataset('Mg_mass',data=Mg_mass)
        grp.create_dataset('Fe_mass',data=Fe_mass)
        grp.create_dataset('MgFe_ratio',data=MgFe_ratio)

    f.close()

if True:
    f = h5py.File(filepath,'r')
    bmap = brewer2mpl.get_map('Spectral','Diverging',9, reverse=True)
    cmap = bmap.mpl_colormap
    color_list = cmap(np.linspace(0.,1.,num=n_snaps))
    # color_list = [self.color_list.append(cmap(0.)),]
    plt.figure()
    for i in np.arange(n_snaps):
        snapnum = snapnum_arr[i]
        print "Plotting snap {}".format(snapnum)

        dat = f['snap{}'.format(str(int(snapnum)).zfill(3))]
        grp_mass = np.array(dat['grp_mass'])
        gas_mass = np.array(dat['gas_mass'])
        Mg_mass = np.array(dat['Mg_mass'])
        Fe_mass = np.array(dat['Fe_mass'])
        MgFe_ratio = np.array(dat['MgFe_ratio'])

        # Set up bins for mass:
        n_mbins = 70
        [mbins_min,mbins_max] = AU._bin_setup(11.,14.,n_mbins)
        mbins_min = 10.**mbins_min
        mbins_max = 10.**mbins_max
        mbins_med = (mbins_min+mbins_max)/2.

        [Q1,med,Q3] = AU._calc_percentiles_v2(grp_mass,mbins_min,mbins_max,MgFe_ratio,min_percentile=20,max_percentile=80)
        good = np.logical_not(np.isnan(Q1))
        plt.plot(np.log10(mbins_min[good]),np.log10(med[good]),color=color_list[i],label='z={}'.format(z_arr[i]))
        plt.fill_between(np.log10(mbins_min[good]),np.log10(Q1[good]),np.log10(Q3[good]),alpha=0.2,color=color_list[i])

    plt.xlabel("Halo Mass")
    plt.ylabel("Log10[Mg/Fe]")
    plt.legend()
    plt.savefig(fig_base+"MgFe_subfind.pdf")





        # plt.savefig(fig_base+"z0_MgFe.pdf")

# save data as npz
