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
parser.add_argument('-tnrm', default=-1, type=int, help='-1 if no timedependent norm, otherwise integer that sets the number of samples which should be in the third column after which to ignore the data')
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
tnrm=args.tnrm

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
tcab,ncab=[],[]
sep2 = sep

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
    ecab = []
    for i in range(sep):
        tcab.append(np.genfromtxt(fcab[i],usecols=1,unpack=True))
        if tnrm != -1: # This part defines what the normalization is. ncab counts if included
            ncab.append((np.genfromtxt(fcab[i],usecols=2,unpack=True)>tnrm)*1)
            tcab[i]=np.multiply(tcab[i],ncab[i])
    energy = np.asarray(energy).astype(float)
    ncab = np.asarray(ncab).astype(float)
    if tnrm != -1:
        sep2 = np.sum(ncab,axis=0)
        for item1 in range(len(inp_n)): # Calcualates the av energy
            tmp=np.multiply(ncab,energy[item1][fstart:fstart+sep][:,None])
            tmpav = np.average(tmp,axis=0)
            ecab.append(tmpav)
    print(len(tcab))
    print(np.shape(np.array(ecab)))
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
            corr = np.divide(corr,sep2,out=np.zeros_like(w1corr),where=sep2!=0)
            w1corr = np.divide(w1corr,sep2,out=np.zeros_like(w1corr),where=sep2!=0)
            w2corr = np.divide(w2corr,sep2,out=np.zeros_like(w1corr),where=sep2!=0)
            w3corr = np.divide(w3corr,sep2,out=np.zeros_like(w1corr),where=sep2!=0)
            w4corr = np.divide(w4corr,sep2,out=np.zeros_like(w1corr),where=sep2!=0)
            """
            corr[:] = [x / float(sep) for x in corr]
            w1corr[:] = [x / float(sep) for x in w1corr]
            w2corr[:] = [x / float(sep) for x in w2corr]
            w3corr[:] = [x / float(sep) for x in w3corr]
            w4corr[:] = [x / float(sep) for x in w4corr]
            """
            # Print Out
            if tnrm != -1:
                np.savetxt('SEG/seg_'+str(int(val))+'_'+item1+'.dat',np.c_[time,ecab[item1count]])
            np.savetxt('SEG/seg_'+str(int(val))+'_'+item1+'_'+item2+'_'+mol_name+'_'+corr_name+'.dat', np.c_[time, corr, w1corr, w2corr, w3corr, w4corr])
            item2count+=1
        item1count+=1

