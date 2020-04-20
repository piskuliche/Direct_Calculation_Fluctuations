#!/usr/bin/env python
import numpy as np
import pickle
import sys,os
import argparse
from scipy import stats
from read_input import user_input

"""
This script has two main roles:
    1) Parse all the correlation functions into something quickly readable
    2) Calculate the correlation function + its derivative.

These two options can be selected with the option flag, and should be run in that order.

Note this script has the support for the first and second derivatives included. Third and fourth derivatives are output as empty files currently, but will be added in at a later date. If these derivatives are needed, then one should use the old programs for calculating these init_segments.py and combine_segments.py.

This script is MUCH faster than its predecessors, especially when considering the ability to read binary files.
It also takes up much less storage space for the same set of calculations.
"""

def calc_avcab(cab, no):
    """
    Calculates the average, with support for a changing normalization
    """
    num = np.sum(cab,axis=0)
    norm = np.sum(no,axis=0)
    return np.divide(num,norm, where=norm!=0, out=np.zeros_like(num))

def zero_order(cab, no):
    """
    Calculates the zeroth order derivative - aka the correlation function!
    Does block averaging, and the error calculation
    """
    d0 = calc_avcab(cab,no)
    nperb = int(nosplit*inputparam.segsplit/inputparam.nblocks)
    bl_d0 = []
    for block in range(inputparam.nblocks):
        bstart=block*nperb
        bend=(block+1)*nperb
        bl_d0.append(calc_avcab(cab[bstart:bend],no[bstart:bend]))
        np.savetxt(subdir+"bl_"+str(block)+"_"+mol_name+'_'+corr_func+'.dat',np.c_[time,bl_d0[block]])
    err_d0 = Error(bl_d0)
    np.savetxt(subdir+item+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,d0,err_d0])
    return

def calc_deriv(cab, encab, no):
    """ 
    Calculates the first derivative
    f' = - < dH f >
    """
    avencab=calc_avcab(encab,no)
    dencab = np.subtract(encab, avencab)
    dcab = -np.multiply(dencab,cab)
    av_dcab = calc_avcab(dcab,nosplit*inputparam.segsplit)
    return av_dcab

def first_order(cab,encab,no):
    """
    This does the first derivative and the blocking!
    """
    d1 = calc_deriv(cab, encab, no)
    nperb = int(nosplit*inputparam.segsplit/inputparam.nblocks)
    bl_d1 = []
    for block in range(inputparam.nblocks):
        bstart=block*nperb
        bend=(block+1)*nperb
        bl_d1.append(calc_deriv(cab[bstart:bend],encab[bstart:bend],no[bstart:bend]))
        np.savetxt(subdir+"bl_"+str(block)+"_"+item+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,bl_d1[block]])
    err_d1 = Error(bl_d1)
    np.savetxt(subdir+item+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,d1,err_d1])
    return 

def calc_deriv2(cab,encab1,encab2,no):
    """
    This calculates the second derivative as:
    f'' = <[dH^2 - <dH^2>]f>
    """
    #Second Derivative
    # <H1>, <H2>
    avencab1, avencab2 = calc_avcab(encab1,no), calc_avcab(encab2,no)
    # dH1 = H1 - <H1>, dH2 = H2-<H2>
    den1cab, den2cab = np.subtract(encab1,avencab1), np.subtract(encab2,avencab2)
    # dH1*dH2 
    d2encab = np.multiply(den1cab,den2cab)
    # <dH1*dH2>
    avd2encab = calc_avcab(d2encab,no)
    # <[dH1*dH2 - <dH1*dH2>]f(t)>
    term1 = calc_avcab(np.multiply(d2encab,cab),nosplit*inputparam.segsplit)
    term2 = calc_avcab(np.multiply(avd2encab,cab),nosplit*inputparam.segsplit)
    av_d2cab = np.subtract(term1,term2)
    return av_d2cab

def second_order(cab,encab,no):
    """
    This does the second derivative for all different options of item
    and blocking and error
    """
    for item2 in inp_n:
        print(item2)
        encab2 = np.multiply(energy[item2][:,None],no)
        d2 = calc_deriv2(cab,encab,encab2,no)
        nperb = int(nosplit*inputparam.segsplit/inputparam.nblocks)
        bl_d2 = []
        for block in range(inputparam.nblocks):
            bstart=block*nperb
            bend=(block+1)*nperb
            bl_d2.append(calc_deriv2(cab[bstart:bend],encab[bstart:bend],encab2[bstart:bend],no[bstart:bend]))
            np.savetxt(subdir+"bl_"+str(block)+"_"+item+'_'+item2+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,bl_d2[block]])
        err_d2 = Error(bl_d2)
        np.savetxt(subdir+item+'_'+item2+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,d2,err_d2])
    return

def third_order(cab,encab,no):
    # Currently not supported.
    d3 = np.zeros(inputparam.corr_length)
    err_d3 = np.zeros(inputparam.corr_length)
    for item2 in inp_n:
        np.savetxt(subdir+item+'_'+item2+'_'+item2+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,d3,err_d3])

def fourth_order(cab,encab,no):
    # Currently not supported.
    d4 = np.zeros(inputparam.corr_length)
    err_d4 = np.zeros(inputparam.corr_length)
    for item2 in inp_n:
        np.savetxt(subdir+item+'_'+item2+'_'+item2+'_'+item2+'_'+mol_name+'_'+corr_func+'.dat',np.c_[time,d4,err_d4])

def Error(bl_corr):
    """
    Calcualtes the uncertainty in the right way from the blocks
    """
    return np.std(bl_corr,axis=0)*t_val



# Reads the command line arguments for the segments
parser = argparse.ArgumentParser()
parser.add_argument('-val', default=0, type=int, help='This is the number of the segment')
parser.add_argument('-opt', default=1, type=int, help='[1] init_segs [2] combine_segs')
parser.add_argument('-corr',default="c2", type=str,help='Correlation function of interest')
parser.add_argument('-mol', default="water",type=str,help="Molecule name")
parser.add_argument('-tnrm',default=-1, type=int, help="[0] regular averaging, [1] special averaging")
parser.add_argument('-fname',default="flucts.inp", type=str, help="File that calls the types of energies")
parser.add_argument('-bin', default=0, type=int, help="[0] Non Binary [1] Binary")
parser.add_argument('-bfcorr', default="norm", type=str, help="norm_ is the default")
parser.add_argument('-time', default="real_time.dat", type=str, help="what file to read the time from")
args = parser.parse_args()
splitno     = args.val
option      = args.opt
corr_func   = args.corr
mol_name    = args.mol
tnrm        = args.tnrm
fname       = args.fname
binaryread  = args.bin
bfcorr      = args.bfcorr
tfile       = args.time

# Read the input file
inputparam = user_input("input_file")

# Pull t-value from Student's T-Table
t_val = stats.t.ppf(0.975,inputparam.nblocks-1)/np.sqrt(inputparam.nblocks)

# Eventually these will be set by read_input
nosplit = int((inputparam.end_config-inputparam.start_config+1000)/(inputparam.segsplit*inputparam.sep_config))

# This calculates the beginning and end of the segment
fstart = splitno*inputparam.sep_config*inputparam.segsplit + inputparam.start_config
fend = (splitno+1)*inputparam.sep_config*inputparam.segsplit + inputparam.start_config

inp_n = np.genfromtxt(fname, usecols=0,dtype=str,unpack=True)
print(inp_n)
subdir="OUT/"


cab,no=[],[]
encab=[]
if option == 1 and tnrm == -1:
    for i in range(fstart,fend,inputparam.sep_config):
        print(i)
        filename="FILES/"+str(i)+"/"+corr_func+"_"+str(i)+"_water"
        # This sets what values above tnrm to contribute, if -1, then disabled
        tmpcab = []
        # Note - support for binary read, much faster!
        if binaryread == 0: tmpcab = np.genfromtxt(filename+".dat", usecols=(1),dtype=(float),unpack=True)
        else: tmpcab = pickle.load(open(filename+".pckl",'rb'))
        # Appends the corr func    
        cab.append(tmpcab)
    if not os.path.exists("TEMP"):
        os.makedirs("TEMP")
    with open('TEMP/cab'+'_'+str(splitno)+'_'+corr_func+'_'+'mol_name'+'.pckl','wb') as g:
        pickle.dump(cab,g)
    # Only if using tnrm
    if tnrm != -1:
        with open('TEMP/no'+'_'+str(splitno)+'_'+corr_func+'_'+'mol_name'+'.pckl','wb') as h:
            pickle.dump(no,h)

elif option == 1 and tnrm != -1:
    for i in range(fstart,fend,inputparam.sep_config):
        print(i)
        filename="FILES/"+str(i)+"/"+corr_func+"_"+str(i)+"_water"
        # This sets what values above tnrm to contribute, if -1, then disabled
        tmpcab,tmpno = [],[]
        # Note - support for binary read - much faster!
        if binaryread == 0: tmpcab,tmpno = np.genfromtxt(filename+'.dat', usecols=(1,2),dtype=(float,float),unpack=True)
        else: 
            fname2 = "FILES/"+str(i)+"/"+bfcorr+"_"+str(i)+"_water"
            tmpcab = pickle.load(open(filename+".pckl",'rb'))
            tmpno  = pickle.load(open(fname2+".pckl", 'rb'))
        no.append((tmpno>tnrm)*1)
        # Appends the corr func
        cab.append(tmpcab)
    if not os.path.exists("TEMP"):
        os.makedirs("TEMP")
    with open('TEMP/cab'+'_'+str(splitno)+'_'+corr_func+'_'+'mol_name'+'.pckl','wb') as g:
        pickle.dump(cab,g)
    # Only if using tnrm
    if tnrm != -1:
        with open('TEMP/no'+'_'+str(splitno)+'_'+corr_func+'_'+'mol_name'+'.pckl','wb') as h:
            pickle.dump(no,h)

elif option == 2:
    # Read in energies and put them in a dictionary
    energy={}
    print("Reading Energies")
    for item in inp_n:
        energy[item]=np.genfromtxt(item+'_init.out',usecols=0)
    # Read in the time, generated by real_time.py
    print("Reading Time")
    time = np.genfromtxt(tfile,usecols=0,unpack=True)
    # This loops over the segments and reads them in. 
    print("Unpickling Files")
    for i in range(nosplit):
        print(i)
        tmpcab = pickle.load(open('TEMP/cab'+'_'+str(i)+'_'+corr_func+'_'+'mol_name'+'.pckl','rb'))
        if tnrm != -1: tmpno = pickle.load(open('TEMP/no'+'_'+str(i)+'_'+corr_func+'_'+'mol_name'+'.pckl','rb'))
        if i == 0:
            cab = tmpcab
            if tnrm != -1: no  = tmpno
        else:
            cab=np.concatenate((cab,tmpcab),axis=0)
            if tnrm != -1: no=np.concatenate((no,tmpno),axis=0)
    # Sets Array for Normalization
    if tnrm == -1: no = np.ones(np.shape(cab))
    # Sets the energy array (note - if tnrm = -1)
    zero_order(cab,no) 
    # Loops over energy items, reads them in to calculate correlations!
    for item in inp_n:
        print("Derivative for  %s" % item)
        encab = np.multiply(energy[item][:,None],no)
        avcab = calc_avcab(cab,no)
        first_order(cab,encab,no)
        second_order(cab,encab,no)
        third_order(cab,encab,no)
        fourth_order(cab,encab,no)
