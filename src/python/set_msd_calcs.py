#set_msd_calcs.py
#This code automatically creates an msd input file for each trajectory folder


from scipy import stats
import numpy as np
import os
import argparse
from argparse import RawTextHelpFormatter
from read_input import input

inputparam = input('../../input_file')
fluctval = int(os.path.split(os.getcwd())[1])


if inputparam.timestep == 'FALSE':
    timestep=1.0
else:
    timestep = float(inputparam.timestep)
with open('mol.info') as f:
        mol_name = str(f.readline().strip())

if inputparam.prog == "LAMMPS":
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
elif inputparam.prog == "CP2K":
    with open("in.nve.cp2k","r") as fi:
        id = []
        for ln in fi:
            if ln.strip().startswith("ABC"):
                id.append(ln.strip().split()[3:])
        vol = id[0][0]


# Writes the new msd calc input file with the right molecule name and volume
filepath='msd_rot_calc.in'
f=open(filepath, 'w')

f.write("# Numeric File\n")
f.write("%s\n" % (fluctval))
f.write("# Number_of_Times Sep_of_Times:\n")
f.write("%s %s\n" % (inputparam.num_times, timestep))
f.write("# Volume\n")
f.write("%s\n" % (vol[0]))
f.write("# Molecule Name\n")
f.write("%s\n" % (mol_name))
f.close()
