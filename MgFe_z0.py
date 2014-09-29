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
CGM_snapbase = "/n/home04/scratch1/AREPOfest/data/CGM_snaps/s120/"

