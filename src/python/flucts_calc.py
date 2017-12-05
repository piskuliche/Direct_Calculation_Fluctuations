# This code is used to calculate the activation energies of diffusion and reorientation based on the fluctuations method.
# It should run for either the NVT or NPT ensembles.

from scipy import stats
import numpy as np
import os
import argparse
from argparse import RawTextHelpFormatter

# Read in Arguments
parser = argparse.ArgumentParser(description='''Calculates the fluctuations of the NPT Ensemble''', formatter_class=RawTextHelpFormatter)
parser.add_argument('-inp', help="Input File")
parser.add_argument('-files', help="Number of Files")
parser.add_argument('-blocks', help="Number of blocks")
args = parser.parse_args()

inputfile = str(args.inp)
nfiles = int(args.files)
nblocks = int(args.blocks)

# Read in input file
inp_names=np.genfromtxt(inputfile, dtype=None,unpack=True)

# Calculate important quantities
filesperblock=nfiles/float(nblocks)
t_val=stats.t.ppf(0.95,nblocks-1)/np.sqrt(nblocks)

# Initialize Vectors
msd=[400]
c2=[400]
ival=[nfiles]
dival=[nfiles]

divalmsd=[400]
divalc1=[400]
divalc2=[400]

eaival_bl=[[] for x in range(0,nblocks)] 
msd_bl=[[] for x in range(0,nblocks)]
dmsdival_bl=[[] for x in range(0,nblocks)]
c2ival_bl=[[] for x in range(0,nblocks)]
c2_bl=[[] for x in range(0,nblocks)]

for i in range(0,nfiles):
    dival.append(0)

# Read in File Names
fnames=np.genfromtxt('file_names', dtype='string', unpack=True)
fmsd=[s+"/msd_"+s+".dat" for s in fnames]
fc1=[s+"/c1_"+s+".dat" for s in fnames]
fc2=[s+"/c2_"+s+".dat" for s in fnames]

# Do the block average calculation
for item in range(0,len(inp_names)):
    print inp_names[item]
    ival=np.genfromtxt(inp_names[item]+'_init.out', unpack=True)
    for block in range(0,nblocks):
        start=int(block*filesperblock)
        end=int((block+1)*filesperblock)
        print start,end
        ivalavg_bl=np.average(ival[start:end])
        for j in range(start,end):
            dival[j]=ival[j]-ivalavg_bl
        for j in range(0, len(msd)):
            msd[j]=0
            divalmsd[j]=0
            c2[j]=0
            divalc2[j]=0
        for j in range(start,end):
            tmsd=np.genfromtxt(fmsd[j], usecols=1, unpack=True)
            tc2=np.genfromtxt(fc2[j], usecols=1, unpack=True)
            # MSD Weighted by Fluct
            msd=msd+tmsd
            divalmsd=divalmsd+tmsd*dival[int(j)]
            # C2 Weighted by Fluct
            c2=c2+tc2
            divalc2=divalc2+tc2*dival[int(j)]
        # Normalize the lists
        # MSD
        msd[:] = [x / float(filesperblock) for x in msd]
        divalmsd[:] = [x / float(filesperblock) for x in divalmsd]
        # C2
        c2[:] = [x / float(filesperblock) for x in c2]
        divalc2[:] = [x / float(filesperblock) for x in divalc2]
        # Calculate Activation Energies
        eaival=divalmsd[1:]/msd[1:]
        # Save to Sep List
        eaival_bl[block] = list(eaival)
        # Saves msd_info to sep list
        msd_bl[block]=list(msd)
        dmsdival_bl[block]=list(divalmsd[1:])
        # Saves c2 to sep list
        c2_bl[block]=list(c2)
        c2ival_bl[block]=list(divalc2)
        # print block values
        np.savetxt('bl_'+str(block)+'_'+inp_names[item]+'_msd.dat', np.c_[eaival], fmt='%s')
        np.savetxt('bl_'+str(block)+'_'+inp_names[item]+'_c2.dat', np.c_[eaival], fmt='%s')
        # Zero the items in the range
        for i in range(start,end):
            dival[i]=0
        for i in range(0,len(msd)):
            msd[i]=0
            divalmsd[i]=0
        for k in range(0,len(eaival)):
            eaival[k]=0
    ivalavg=np.average(ival)
    print ivalavg
    # Zeroes the msd vectors and recalculates the d_values
    for j in range(0,nfiles):
        dival[j]=ival[j]-ivalavg
    for j in range(0,len(msd)):
        msd[j]=0
        divalmsd[j]=0
        c2[j]=0
        divalc2[j]=0
    for j in range(0,nfiles):
        tmsd=np.genfromtxt(fmsd[j], usecols=1, unpack=True)
        tc2=np.genfromtxt(fc2[j], usecols=1, unpack=True)
        # Weighted By Flucts
        # MSD
        msd=msd+tmsd
        divalmsd=divalmsd+tmsd*dival[int(j)]
        # C2
        c2=c2+tc2
        divalc2=divalc2+tc2*dival[int(j)]
    # Normalize Results
    msd[:] = [x / float(nfiles) for x in msd]
    divalmsd[:] = [x / float(nfiles) for x in divalmsd]
    c2[:] = [x / float(nfiles) for x in c2]
    divalc2 = [x / float(nfiles) for x in divalc2]
    eaival=divalmsd[1:]/msd[1:]
    # Save Data
    np.savetxt('d'+inp_names[item]+'.dat', dival, fmt=['%.4f'])
    # Calculate Uncertainty
    # MSD
    err_ivalmsdea=np.array(eaival_bl).std(0)
    err_ivalmsdea= [x * t_val for x in err_ivalmsdea]
    # DMSD
    err_divalmsd=np.array(dmsdival_bl).std(0)
    err_divalmsd= [x * t_val for x in err_divalmsd]
    # DC2
    err_c2ival=np.array(c2ival_bl).std(0)
    err_c2ival= [x * t_val for x in err_c2ival]

    # Calculate Time
    time = []
    for i in range(0, len(err_ivalmsdea)+1):
        time.append(i*0.05)
    print len(eaival)
    print len(time)
    print len(err_ivalmsdea)
    print len(divalc2)
    print len(err_c2ival)
    print len(c2)

    np.savetxt('ea_msd_'+inp_names[item]+'.dat', np.c_[time[1:], eaival, err_ivalmsdea], fmt='%s')
    np.savetxt('d'+inp_names[item]+'_msd_tot_.dat', np.c_[time[1:], divalmsd[1:], err_divalmsd], fmt='%s')
    np.savetxt('dc2_'+inp_names[item]+'.dat', np.c_[time, divalc2, err_c2ival], fmt='%s')
    np.savetxt('c2_total_result.dat', np.c_[time,c2], fmt='%s')
    np.savetxt('dc2_'+inp_names[item]+'_total_result.dat', np.c_[time, divalc2], fmt='%s')
    np.savetxt('time.dat',time,fmt='%s')
    for i in range(0, nblocks):
        np.savetxt('bl_'+str(int(i))+'_c2_val.dat', c2_bl[i], fmt='%s')
        np.savetxt('bl_'+str(int(i))+'_dc2_val_'+inp_names[item]+'.dat', c2ival_bl[i], fmt='%s')
