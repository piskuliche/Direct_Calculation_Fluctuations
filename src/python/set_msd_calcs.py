#set_msd_calcs.py
#This code automatically creates an msd input file for each trajectory folder


from scipy import stats
import numpy as np
import os
import argparse
from argparse import RawTextHelpFormatter

# Read Input from Command Line
parser = argparse.ArgumentParser(description='''Sets up the msd calculations''', formatter_class=RawTextHelpFormatter)
parser.add_argument('-inp', help="Input File Number")
parser.add_argument('-stp', help="Timestep (in ps)")
parser.add_argument('-ntimes', help="Number of distinct times")
parser.add_argument('-mol', help="Molecule Name")
args = parser.parse_args()
fluctval = int(args.inp)
if args.stp == 'FALSE':
    timestep=1.0
else:
    timestep = float(args.stp)

ntimes   = int(args.ntimes)
mol_name = str(args.mol)

# Looks up the Step and Loop time keywords from log.lammps
filename='log.lammps'
lookup= 'Step'
lookup2= 'Loop time'
totlines=0
endskip=0
with open(filename) as myFile:
    for num, line in enumerate(myFile, 1):
        if lookup in line:
            startskip=num
        if lookup2 in line:
            endskip=num
        totlines=num
endskip=totlines-endskip-4
vol = np.genfromtxt(filename, skip_header=startskip,skip_footer=endskip,usecols=(11), unpack=True)

# Writes the new msd calc input file with the right molecule name and volume
filepath='msd_rot_calc.in'
f=open(filepath, 'w')
f.write("# Numeric File\n")
f.write("%s\n" % (fluctval))
f.write("# Number_of_Times Sep_of_Times:\n")
f.write("%s %s\n" % (ntimes, timestep))
f.write("# Volume\n")
f.write("%s\n" % (vol[0]))
f.write("# Molecule Name\n")
f.write("%s\n" % (mol_name))
f.close()
