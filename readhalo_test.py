import readhaloHDF5

snapnum=120
snapbase="/n/ghernquist/Illustris/Runs/Illustris-1/output/"
type=0
grpnr=80
subnr=-1
readhaloHDF5.reset()
mass=readhaloHDF5.readhalo(snapbase, "snap", snapnum, "MASS", type, grpnr, subnr, long_ids=True, double_output=False)

print "mass ",mass