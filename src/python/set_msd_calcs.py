from scipy import stats
import numpy as np
import os
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description='''Sets up the msd calculations''', formatter_class=RawTextHelpFormatter)
parser.add_argument('-inp', help="Input File")
args = parser.parse_args()

fluctval = int(args.inp)

filename='log.lammps'
lookup= 'Step'
lookup2= 'Loop time'
totlines=0
with open(filename) as myFile:
    for num, line in enumerate(myFile, 1):
        if lookup in line:
            startskip=num
        if lookup2 in line:
            endskip=num
        totlines=num
endskip=totlines-endskip-4
vol = np.genfromtxt(filename, skip_header=startskip,skip_footer=endskip,usecols=(11), unpack=True)

filepath='msd_rot_calc.in'
f=open(filepath, 'w')
f.write("# Numeric File\n")
f.write("%s\n" % (fluctval))
f.write("# Number_of_Times Sep_of_Times:\n")
f.write("400 0.050\n")
f.write("# Volume\n")
f.write("%s\n" % (vol[0]))
f.close()
