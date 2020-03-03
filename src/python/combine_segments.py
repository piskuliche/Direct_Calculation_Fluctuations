#!/usr/bin/env python
import numpy as np
import sys
import argparse
from scipy import stats
from read_input import user_input


"""
This is a python program to take the weighted correlation functions calculated in init_segments.py and turn them into the derivative correlation functions.
"""

def BLOCK_ENERGY(energy,start,end,itemindex):
    eav = 0
    for i in range(start*sep, end*sep):
        eav += energy[itemindex][i]
    eav = eav/(end*sep-start*sep)
    return eav

def BLOCK_DENERGY(energy,start,end,itemindex,eav,n):
    deav = 0
    for i in range(start*sep, end*sep):
        deav += (energy[itemindex][i]-eav)**n
    deav = deav/(end*sep-start*sep)
    return deav

def BLOCK_ARRAY_ERR(barray, noblocks):
    sigma = []
    for i in range(len(barray)):
        sigma.append(np.std(barray[i])*t_val)
    return sigma

def NORM(array, n):
    return array/float(n)

"""
Subtracts the averages:
    Follows identities of deltaH from pascal's triangle
    dH(0) = H-<H>
    dH(0)^2 = H^2 - 2*H*<H> + <H>^2
    dH(0)^3 = H^3 - 3*H^2*<H> + 3*H*<H>^2 - <H>^3
    dH(0)^4 = H^4 - 4*H^3*<H> + 6*H^2*<H>^2 - 4*H*<H>^3+<H>^4
"""
def FRST_SUB_AV(corr,w1corr,av):
    w1corr = w1corr-av*corr
    return w1corr

def SCND_SUB_AV(corr,w1corr,w2corr,av):
    w2corr = w2corr - 2*w1corr*av + (av**2.)*corr
    return w2corr

def THRD_SUB_AV(corr,w1corr,w2corr,w3corr,av):
    w3corr = w3corr - 3*w2corr*av + 3*w1corr*(av**2.) - (av**3.)*corr
    return w3corr

def FRTH_SUB_AV(corr,w1corr,w2corr,w3corr,w4corr,av):
    w4corr = w4corr - 4*w3corr*av + 6*w2corr*(av**2.) - 4*w1corr*(av**3) + (av**4.)*corr
    return w4corr

    
"""
Derivative for the first derivative
C'(t) = -<dh(0)A(0)B(t)>
"""
def FIRST_DERIV(w1corr):
    d1corr = -w1corr
    return d1corr

"""
Derivative for the second derivative
C''(t) = <dh(0)^2A(0)B(t)>-<dh^2>C(t)
"""
def SECOND_DERIV(corr, w2corr, d2av):
    d2corr = w2corr - d2av*corr
    return d2corr

"""
Derivative for the third derivative
C'''(t) = -<dh(0)^3A(0)B(t)> + <dh^3>C(t) - 3<dh^2>C'(t)
"""
def THIRD_DERIV(corr, d1corr, w3corr, d2av, d3av):
    d3corr = -w3corr - 3*d2av*d1corr + d3av*corr
    return d3corr
"""
Derivative for the fourth derivative
C''''(t) = <[dh(0)^4 -<dh^4>]A(0)B(t)> - 6*<dh^2>C''(t)+4<dh^3>C'(t)
"""
def FOURTH_DERIV(corr, d1corr, d2corr, w4corr, d2av, d3av, d4av):
    d4corr = w4corr - 6*d2av*d2corr + 4*d3av*d1corr - d4av*corr
    return d4corr

def RATIO(corr, d1corr):
    ea = np.divide(d1corr[1:],corr[1:])
    ea=np.insert(ea,0,0.0)
    return ea
    


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

# Initialize Arrays
time   = np.zeros(corrlength)
blcorr   = np.zeros((inputparam.nblocks,len(inp_n),len(inp_n), num_segs, corrlength))
blw1corr = np.zeros((inputparam.nblocks,len(inp_n),len(inp_n), num_segs, corrlength))
blw2corr = np.zeros((inputparam.nblocks,len(inp_n),len(inp_n), num_segs, corrlength))
blw3corr = np.zeros((inputparam.nblocks,len(inp_n),len(inp_n), num_segs, corrlength))
blw4corr = np.zeros((inputparam.nblocks,len(inp_n),len(inp_n), num_segs, corrlength))
corr   = np.zeros((len(inp_n),len(inp_n), num_segs, corrlength))
w1corr = np.zeros((len(inp_n),len(inp_n), num_segs, corrlength))
w2corr = np.zeros((len(inp_n),len(inp_n), num_segs, corrlength))
w3corr = np.zeros((len(inp_n),len(inp_n), num_segs, corrlength))
w4corr = np.zeros((len(inp_n),len(inp_n), num_segs, corrlength))

energy = []

# Read Segments In
item1count=0
print("---Start Read---")
for item1 in inp_n:
    print("The current Item is: %s" % item1)
    sys.stdout.flush()
    item2count=0
    for item2 in inp_n:
        for seg in range(num_segs):
            # 4D Arrays of shape [item1][item2][seg][time]
            time,corr[item1count][item2count][seg],w1corr[item1count][item2count][seg],w2corr[item1count][item2count][seg],w3corr[item1count][item2count][seg],w4corr[item1count][item2count][seg]=np.genfromtxt('SEG/seg_'+str(seg)+'_'+item1+'_'+item2+'_'+mol_name+'_'+corr_name+'.dat', usecols=(0,1,2,3,4,5), unpack=True)
        item2count += 1
    item1count += 1
print("---End Read---")
print("---Start Blocking---")
for item1 in inp_n:
    energy.append(np.genfromtxt(item1+'_init.out'))
item1count = 0
for item1 in inp_n:
    # This overrides the average calc and reads it from file
    if inputparam.cab=="IONPAIRING":
        with open(corr_name+'_react_'+str(item1)+'.dat','r') as f:
            reav=float(f.readline())
    print("Starting loop for %s" % item1)
    item2count = 0
    for item2 in inp_n:
        print("Calculating %s %s" % (item1,item2))
        item3 = item2
        item4 = item2
        # Zero Block Weighted Arrays
        bl_w1corr  = np.zeros((inputparam.nblocks, corrlength))
        bl_w2corr  = np.zeros((inputparam.nblocks, corrlength))
        bl_w3corr  = np.zeros((inputparam.nblocks, corrlength))
        bl_w4corr  = np.zeros((inputparam.nblocks, corrlength))
        # Zero Block Derivative Arrays 
        bl_corr    = np.zeros((inputparam.nblocks, corrlength))
        bl_d1corr  = np.zeros((inputparam.nblocks, corrlength))
        bl_d2corr  = np.zeros((inputparam.nblocks, corrlength))
        bl_d3corr  = np.zeros((inputparam.nblocks, corrlength))
        bl_d4corr  = np.zeros((inputparam.nblocks, corrlength))
        bl_ea      = np.zeros((inputparam.nblocks, corrlength))
        # Zero Tot Weighted Arrays
        tot_w1corr = np.zeros(corrlength)
        tot_w2corr = np.zeros(corrlength)
        tot_w3corr = np.zeros(corrlength)
        tot_w4corr = np.zeros(corrlength)
        # Zero Tot Derivative Arrays
        tot_corr   = np.zeros(corrlength)
        tot_d1corr = np.zeros(corrlength)
        tot_d2corr = np.zeros(corrlength)
        tot_d3corr = np.zeros(corrlength)
        tot_d4corr = np.zeros(corrlength)
        tot_ea     = np.zeros(corrlength)
        # Zero Err Weighted Arrays
        err_w1corr = np.zeros(corrlength)
        err_w2corr = np.zeros(corrlength)
        err_w3corr = np.zeros(corrlength)
        err_w4corr = np.zeros(corrlength)
        # Zero Err Derivative Arrays
        err_corr   = np.zeros(corrlength)
        err_d1corr = np.zeros(corrlength)
        err_d2corr = np.zeros(corrlength)
        err_d3corr = np.zeros(corrlength)
        err_d4corr = np.zeros(corrlength)
        err_ea     = np.zeros(corrlength)
        # Sum Segments into blocks
        for block in range(inputparam.nblocks):
            print("     BLOCK %s" % block)
            # Block indices calculation
            bstart = int(block*segs_per_block)
            bend   = int((block+1)*segs_per_block)
            bdist  = bend - bstart
            # Calculate Average Flucts
            if inputparam.cab=="TRANSPORT":
                e1_av  =   BLOCK_ENERGY(energy, bstart, bend, item1count)
            elif inputparam.cab=="IONPAIRING":
                e1_av = reav
            else:
                print("Error: Type not TRANSPORT or IONPAIRING")
                exit()
            d1_av  =  BLOCK_DENERGY(energy, bstart, bend, item1count,e1_av,1)
            d2_av  =  BLOCK_DENERGY(energy, bstart, bend, item1count,e1_av,2)
            d3_av  =  BLOCK_DENERGY(energy, bstart, bend, item1count,e1_av,3)
            d4_av  =  BLOCK_DENERGY(energy, bstart, bend, item1count,e1_av,4)
            print("-------")
            print("e1av = %s" % e1_av)
            print("d1av = %s" % d1_av)
            print("d2av = %s" % d2_av)
            print("d3av = %s" % d3_av)
            print("d4av = %s" % d4_av)
            for i in range(corrlength):
                for seg in range(bstart,bend):
                    # Need to normalize.
                    bl_corr[block][i]   +=   corr[item1count][item2count][seg][i]
                    bl_w1corr[block][i] += FRST_SUB_AV(corr[item1count][item2count][seg][i], w1corr[item1count][item2count][seg][i],e1_av)
                    bl_w2corr[block][i] += SCND_SUB_AV(corr[item1count][item2count][seg][i], w1corr[item1count][item2count][seg][i],w2corr[item1count][item2count][seg][i],e1_av)
                    bl_w3corr[block][i] += THRD_SUB_AV(corr[item1count][item2count][seg][i], w1corr[item1count][item2count][seg][i],w2corr[item1count][item2count][seg][i],w3corr[item1count][item2count][seg][i],e1_av)
                    bl_w4corr[block][i] += FRTH_SUB_AV(corr[item1count][item2count][seg][i], w1corr[item1count][item2count][seg][i],w2corr[item1count][item2count][seg][i],w3corr[item1count][item2count][seg][i],w4corr[item1count][item2count][seg][i],e1_av)
                # Average over segments
                bl_corr[block][i]   = NORM(bl_corr[block][i], bdist)
                bl_w1corr[block][i] = NORM(bl_w1corr[block][i], bdist)
                bl_w2corr[block][i] = NORM(bl_w2corr[block][i], bdist)
                bl_w3corr[block][i] = NORM(bl_w3corr[block][i], bdist)
                bl_w4corr[block][i] = NORM(bl_w4corr[block][i], bdist)
                # Calculate Derivatives
                bl_d1corr[block][i] = FIRST_DERIV(bl_w1corr[block][i])
                bl_d2corr[block][i] = SECOND_DERIV(bl_corr[block][i],bl_w2corr[block][i],d2_av)
                bl_d3corr[block][i] = THIRD_DERIV(bl_corr[block][i], bl_d1corr[block][i],bl_w3corr[block][i],d2_av, d3_av)
                bl_d4corr[block][i] = FOURTH_DERIV(bl_corr[block][i], bl_d1corr[block][i], bl_d2corr[block][i], bl_w4corr[block][i], d2_av, d3_av, d4_av)
            # Calculate ratio function
            bl_ea[block] = RATIO(bl_corr[block], bl_d1corr[block])
            # Sets File Names
            bl_name   = "bl_"+str(block)+"_"+mol_name+"_"+corr_name+".dat"
            bl_d1name = "bl_"+str(block)+"_"+item1+"_"+mol_name+"_"+corr_name+".dat"
            bl_d2name = "bl_"+str(block)+"_"+item1+"_"+item2+"_"+mol_name+"_"+corr_name+".dat"
            bl_d3name = "bl_"+str(block)+"_"+item1+"_"+item2+"_"+item3+"_"+mol_name+'_'+corr_name+".dat"
            bl_d4name = "bl_"+str(block)+"_"+item1+"_"+item2+"_"+item3+"_"+item4+"_"+mol_name+'_'+corr_name+".dat"
            # Writes to file
            np.savetxt(bl_name,   np.c_[time, bl_corr[block], bl_ea[block]],      fmt='%s')
            np.savetxt(bl_d1name, np.c_[time, bl_d1corr[block],bl_w1corr[block]], fmt='%s')
            np.savetxt(bl_d2name, np.c_[time, bl_d2corr[block],bl_w2corr[block]], fmt='%s')
            np.savetxt(bl_d3name, np.c_[time, bl_d3corr[block],bl_w3corr[block]], fmt='%s')
            np.savetxt(bl_d4name, np.c_[time, bl_d4corr[block],bl_w4corr[block]], fmt='%s')
        # Calculate Uncertainties
        err_corr    =   np.array(bl_corr).std(0)
        err_corr    =   [x * t_val for x in err_corr]
        err_w1corr  =   np.array(bl_w1corr).std(0)
        err_w1corr  =   [x * t_val for x in err_w1corr]
        err_w2corr  =   np.array(bl_w2corr).std(0)
        err_w2corr  =   [x * t_val for x in err_w2corr]
        err_w3corr  =   np.array(bl_w3corr).std(0)
        err_w3corr  =   [x * t_val for x in err_w3corr]
        err_w4corr  =   np.array(bl_w4corr).std(0)
        err_w4corr  =   [x * t_val for x in err_w4corr]
        err_d1corr  =   np.array(bl_d1corr).std(0)
        err_d1corr  =   [x * t_val for x in err_d1corr]
        err_d2corr  =   np.array(bl_d2corr).std(0)
        err_d2corr  =   [x * t_val for x in err_d2corr]
        err_d3corr  =   np.array(bl_d3corr).std(0)
        err_d3corr  =   [x * t_val for x in err_d3corr]
        err_d4corr  =   np.array(bl_d4corr).std(0)
        err_d4corr  =   [x * t_val for x in err_d4corr]
        err_ea      =   np.array(bl_ea).std(0)
        err_ea      =   [x * t_val for x in err_corr]
        # Calculate Averages
        seg_start = 0
        seg_end   = num_segs
        seg_dist  = seg_end - seg_start
        if inputparam.cab=="TRANSPORT":
            e1_av  =   BLOCK_ENERGY(energy, seg_start, seg_end, item1count)
        elif inputparam.cab=="IONPAIRING":
            e1_av = reav
        else:
            print("Error: Type not TRANSPORT or IONPAIRING")
            exit()
        e1_av = BLOCK_ENERGY(energy, seg_start, seg_end, item1count)
        d1_av = BLOCK_DENERGY(energy, seg_start, seg_end, item1count,e1_av, 1)
        d2_av = BLOCK_DENERGY(energy, seg_start, seg_end, item1count,e1_av, 2)
        d3_av = BLOCK_DENERGY(energy, seg_start, seg_end, item1count,e1_av, 3)
        d4_av = BLOCK_DENERGY(energy, seg_start, seg_end, item1count,e1_av, 4)
        print("e1_av = %s" % e1_av)
        print("d1_av = %s" % d1_av)
        print("d2_av = %s" % d2_av)
        print("d3_av = %s" % d3_av)
        print("d4_av = %s" % d4_av)
        # Loop over number of times
        for i in range(corrlength):
            for seg in range(num_segs):
                tot_corr[i]   += corr[item1count][item2count][seg][i]
                tot_w1corr[i] += FRST_SUB_AV(corr[item1count][item2count][seg][i], w1corr[item1count][item2count][seg][i],e1_av)
                tot_w2corr[i] += SCND_SUB_AV(corr[item1count][item2count][seg][i], w1corr[item1count][item2count][seg][i],w2corr[item1count][item2count][seg][i],e1_av)
                tot_w3corr[i] += THRD_SUB_AV(corr[item1count][item2count][seg][i], w1corr[item1count][item2count][seg][i],w2corr[item1count][item2count][seg][i],w3corr[item1count][item2count][seg][i],e1_av)
                tot_w4corr[i] += FRTH_SUB_AV(corr[item1count][item2count][seg][i], w1corr[item1count][item2count][seg][i],w2corr[item1count][item2count][seg][i],w3corr[item1count][item2count][seg][i],w4corr[item1count][item2count][seg][i],e1_av)
            # Normalize
            tot_corr[i] = NORM(tot_corr[i], seg_dist) 
            tot_w1corr[i] = NORM(tot_w1corr[i], seg_dist)
            tot_w2corr[i] = NORM(tot_w2corr[i], seg_dist)
            tot_w3corr[i] = NORM(tot_w3corr[i], seg_dist)
            tot_w4corr[i] = NORM(tot_w4corr[i], seg_dist)
            # Calculate Derivatives
            tot_d1corr[i] = FIRST_DERIV(tot_w1corr[i])
            tot_d2corr[i] = SECOND_DERIV(tot_corr[i],tot_w2corr[i],d2_av)
            tot_d3corr[i] = THIRD_DERIV(tot_corr[i], tot_d1corr[i],tot_w3corr[i],d2_av, d3_av)
            tot_d4corr[i] = FOURTH_DERIV(tot_corr[i], tot_d1corr[i], tot_d2corr[i], tot_w4corr[i], d2_av, d3_av, d4_av)
        # Calculate ratio of correlation functions
        tot_ea = RATIO(tot_corr,tot_d1corr)
        # Name Correlation Functions
        tot_name   = mol_name+"_"+corr_name+".dat"
        tot_d1name = item1+"_"+mol_name+"_"+corr_name+".dat"
        tot_d2name = item1+"_"+item2+"_"+mol_name+"_"+corr_name+".dat"
        tot_d3name = item1+"_"+item2+"_"+item3+"_"+mol_name+"_"+corr_name+".dat"
        tot_d4name = item1+"_"+item2+"_"+item3+"_"+item4+"_"+mol_name+"_"+corr_name+".dat"
        # Print Correlation Functions
        np.savetxt(tot_name,   np.c_[time, tot_corr, err_corr, tot_ea, err_ea], fmt='%s')
        np.savetxt(tot_d1name, np.c_[time, tot_d1corr, err_d1corr, tot_w1corr, err_w1corr], fmt='%s')
        np.savetxt(tot_d2name, np.c_[time, tot_d2corr, err_d2corr, tot_w2corr, err_w2corr], fmt='%s')
        np.savetxt(tot_d3name, np.c_[time, tot_d3corr, err_d3corr, tot_w3corr, err_w3corr], fmt='%s')
        np.savetxt(tot_d4name, np.c_[time, tot_d4corr, err_d4corr, tot_w4corr, err_w4corr], fmt='%s')
        
        item2count+=1
    item1count+=1


            
            




               







