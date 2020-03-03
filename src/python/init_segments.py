#!/usr/bin/env python
import numpy as np
import sys
import argparse
from read_input import user_input

def BLOCK_ENERGY(energy,start,end,itemindex):
    eav = 0
    for i in range(start*sep, end*sep):
        eav += energy[itemindex][i]
    eav = eav/(end*sep-start*sep)
    return eav

# Override Parameters
parser = argparse.ArgumentParser()
parser.add_argument('-val', default=0, type=int, help='This is the number of the segment')
parser.add_argument('-fname', default="flucts.inp", type=str, help='This sets the fluctuations file name')
parser.add_argument('-corr', default="c2", type=str, help='This is the correlation function name')
parser.add_argument('-mol', default="water", type=str, help='This is the molecule name')
parser.add_argument('-timeoverride', default=0, type=int, help='This overrides the time to be something other than the time')
parser.add_argument('-foverride', default="time.override", type=str, help='This is the file the time override is read from')
args = parser.parse_args()
usetime = args.timeoverride
foverride = args.foverride
val = args.val
fname=args.fname
corr_name=args.corr
mol_name = args.mol

# Read the input file
inputparam = user_input("input_file")

# Set the number of segments
sep = int(inputparam.segsplit)
num_segs = int(inputparam.num_files/float(sep))
segs_per_block = np.floor(num_segs/float(inputparam.nblocks))


# Read In Time
time = []
corrlength = 0
if usetime == 0:
    time = np.genfromtxt('real_time.dat', usecols=0)
    time = [0] + time
    time = time[:inputparam.num_times]
    corrlength = inputparam.num_times
else:
    time = np.genfromtxt(foverride, usecols=0)
    corrlength = len(time)

# Read In Fluctuations
inp_n = np.genfromtxt(fname, usecols=0, dtype=str,unpack=True)

# Initialize Arrays
d1=[sep]
d2=[sep]
d3=[sep]
d4=[sep]
energy=[]
tcab=[]

print("The Argument Provided is %s:" % val)
if val == '-h':
    print("This function needs a single integer input to run.")
    print("     Usage: python init_average.py sep fname corr_func")
    exit()
else:
    itmp=int(val)
    print(itmp)
    fstart = itmp*sep
    fend = (itmp+1)*sep
    print("Starting Trajectory: %s\nEnding Trajectory: %s" % (fstart,fend))
    # Read in Correlation Functions
    fnames = np.genfromtxt('file_names', dtype=str, unpack=True)
    print("Choosing %s of %s files" % (sep, len(fnames)))
    # Chooses where to pull from in cases when there are a lot of files.
    if itmp < 500:
        fcab=["FILES/"+str(s)+"/"+str(corr_name)+"_"+str(s)+"_"+str(mol_name)+".dat" for s in fnames[fstart:fend]]
    else:
        fcab=["FILES2/"+str(s)+"/"+str(corr_name)+"_"+str(s)+"_"+str(mol_name)+".dat" for s in fnames[fstart:fend]]
    item1count=0
    # Read in the energy files
    for item1 in inp_n:
        energy.append(np.genfromtxt(item1+'_init.out'))
        print("Average of %s is: %s" % (item1, np.average(energy[item1count])))
        item1count += 1
    item1count = 0
    for i in range(sep):
        tcab.append(np.genfromtxt(fcab[i],usecols=1,unpack=True))
    print(len(tcab))
    for item1 in inp_n:
        print(item1)
        item2count=0
        for item2 in inp_n:
            # Calculate Average Energy
            seg_start = 0
            seg_end = num_segs
            seg_dist = seg_end - seg_start
            # Initialize correlation functions
            corr=np.zeros((corrlength))
            w1corr=np.zeros((corrlength))
            w2corr=np.zeros((corrlength))
            w3corr=np.zeros((corrlength))
            w4corr=np.zeros((corrlength))
            # Loop over files in segment
            for i in range(sep):
                enum=fstart+i
                corr=corr+tcab[i]
                d1=energy[item1count][enum]
                d2=energy[item2count][enum]
                w1corr=w1corr+tcab[i]*d1
                w2corr=w2corr+tcab[i]*d1*d2
                w3corr=w3corr+tcab[i]*d1*d2*d2
                w4corr=w4corr+tcab[i]*d1*d2*d2*d2
            # Normalize
            corr[:] = [x / float(sep) for x in corr]
            w1corr[:] = [x / float(sep) for x in w1corr]
            w2corr[:] = [x / float(sep) for x in w2corr]
            w3corr[:] = [x / float(sep) for x in w3corr]
            w4corr[:] = [x / float(sep) for x in w4corr]
            # Print Out
            np.savetxt('SEG/seg_'+str(int(val))+'_'+item1+'_'+item2+'_'+mol_name+'_'+corr_name+'.dat', np.c_[time, corr, w1corr, w2corr, w3corr, w4corr])
            item2count+=1
        item1count+=1

