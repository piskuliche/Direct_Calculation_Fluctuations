#!/usr/bin/env python
import numpy as np
import sys
import argparse
from scipy import stats
from read_input import user_input

# Import Input File Parameters
inputparam = user_input('input_file')

# Override Parameters
parser = argparse.ArgumentParser()
parser.add_argument('-fname', default="flucts.inp", type=str, help='This sets the fluctuations file name')
parser.add_argument('-corr', default="c2", type=str, help='This is the correlation function name')
parser.add_argument('-mol', default="water", type=str, help='This is the molecule name')
parser.add_argument('-timeoverride', default=0, type=int, help='This overrides the time to be something other than the time')
parser.add_argument('-foverride', default="time.override", type=str, help='This is the file the time override is read from')
args = parser.parse_args()
usetime = args.timeoverride
foverride = args.foverride
fname = args.fname
corr_name = args.corr
mol_name = args.mol

corrlength=0
if usetime == 0:
    corrlength = inputparam.num_times
else:
    t = np.genfromtxt(foverride, usecols=0)
    corrlength = len(t)

# Pull t-value from Student's T-Table
t_val = stats.t.ppf(0.975,inputparam.nblocks-1)/np.sqrt(inputparam.nblocks)

#Calculate number of segments
sep = int(inputparam.segsplit)
num_segs = int(inputparam.num_files/float(sep))
segs_per_block = np.floor(num_segs/float(inputparam.nblocks))
print("There are %s total segments" % num_segs)
print("There are %s blocks" % inputparam.nblocks)
print("There are %s segs_per_block" % segs_per_block)

# Read In Fluctuations
inp_n = np.genfromtxt(fname, usecols=0, dtype=str,unpack=True)

seg_en = np.zeros((len(inp_n), num_segs, corrlength))
seg_no = np.zeros((len(inp_n), num_segs, corrlength))
corr = np.zeros((len(inp_n), num_segs, corrlength))
time  = np.zeros(corrlength)
wcorr = np.zeros((len(inp_n), num_segs, corrlength))
dcorr = np.zeros((len(inp_n), num_segs, corrlength))
b_enav = np.zeros((inputparam.nblocks,len(inp_n),corrlength))
b_corr = np.zeros((inputparam.nblocks,len(inp_n),corrlength))
b_wcorr = np.zeros((inputparam.nblocks,len(inp_n),corrlength))
b_dcorr = np.zeros((inputparam.nblocks,len(inp_n),corrlength))

index = 0
print("Starting Read")
for item in inp_n:
    for seg in range(num_segs):
        seg_en[index][seg] = np.genfromtxt("SEG/seg_"+str(seg)+'_'+item+'.dat',usecols=1)
        seg_no[index][seg] = (seg_en[index][seg]!=0)*1
        time, corr[index][seg], wcorr[index][seg] = np.genfromtxt("SEG/seg_"+str(seg)+"_"+item+"_"+item+"_"+mol_name+'_'+corr_name+".dat",usecols=(0,1,2),unpack=True)
    index +=1
print("Read is Complete")
print("Starting Blocking")
"""
for block in range(inputparam.nblocks):
    print("Block %d" % block)
    # Block indices calculation
    bstart = int(block*segs_per_block)
    bend   = int((block+1)*segs_per_block)
    bdist  = bend - bstart
    print(np.shape(seg_en),np.shape(seg_no))
    num_en, div_en = np.sum(seg_en[:,bstart:bend,:],axis=1), np.sum(seg_no[:,bstart:bend,:],axis=1)
    print(np.shape(num_en),np.shape(div_en))
    b_enav[block]=np.divide(num_en,div_en, where=div_en!=0, out=np.zeros_like(num_en))
    num_corr, div_corr = np.sum(corr[:,bstart:bend,:],axis=1), np.sum(corr[:,bstart:bend,:],axis=1)
    b_corr[block]=np.divide(num_corr,div_corr, where=div_corr!=0, out=np.zeros_like(num_corr))
    num_wcorr, div_wcorr = np.sum(wcorr[:,bstart:bend,:],axis=1), np.sum(wcorr[:,bstart:bend,:],axis=1)
    b_wcorr[block]=np.divide(num_corr,div_corr, where=div_corr!=0, out=np.zeros_like(num_corr))
    b_dcorr[block]=np.subtract(b_wcorr[block],np.multiply(b_corr[block],b_enav[block]))
print("Blocking done")
"""
num_en, div = np.sum(seg_en,axis=1), np.sum(seg_no,axis=1)
enav=np.divide(num_en,div, where=div!=0, out=np.zeros_like(num_en))
num_corr = np.sum(corr,axis=1)
fcorr=np.divide(num_corr,div, where=div!=0, out=np.zeros_like(num_corr))
term2=np.multiply(corr,enav[:,None,:])
fwcorr=-np.subtract(wcorr,term2)
num_wcorr = np.sum(fwcorr,axis=1)
fdcorr=np.divide(num_wcorr,div, where=div!=0, out=np.zeros_like(num_wcorr))

np.savetxt('testout',np.c_[time,fcorr[0],fdcorr[0],enav[0]])


