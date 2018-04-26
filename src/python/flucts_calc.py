#!/usr/bin/env python
# This code is used to calculate the activation energies of diffusion and reorientation based on the fluctuations method.
# It should run for either the NVT or NPT ensembles.

from scipy import stats
import numpy as np
import os
import argparse
from argparse import RawTextHelpFormatter
import sys


def Third_Order_Weight(i, d3, d3av):
    return -(d3[i]-d3av)

def Fourth_Order_Weight(i,d4,d4av):
    return (d4[i] - d4av)


# Read in Arguments
parser = argparse.ArgumentParser(description='''Calculates the fluctuations of the NPT Ensemble''', formatter_class=RawTextHelpFormatter)
parser.add_argument('-inp', help="Input File")
parser.add_argument('-files', help="Number of Files")
parser.add_argument('-blocks', help="Number of blocks")
parser.add_argument('-ntimes', help="Number of times")
parser.add_argument('-mol', help ="Molecule Name")
parser.add_argument('-ind', help = "If you want a single element - type its position in the input file (counting from 0). If you want the whole thing - type FALSE")
args = parser.parse_args()

inputfile = str(args.inp) # This is the input file that specifies the items and the columns
nfiles = int(args.files) # This is the total number of files
nblocks = int(args.blocks) # This is the number of blocks for averaging
ntimes = int(args.ntimes) # This is the length of the TCF calculation
mol_name = str(args.mol) # This is the name of the molecule


inp_n, inp_c=np.genfromtxt(inputfile, usecols=(0,1), dtype=(str,int),unpack=True)
index=0
inp_names1=[]
inp_names2=[]
print args.ind
if args.ind == 'FALSE' or args.ind == 'None':
    inp_names1=inp_n
    inp_names2=inp_n
    inp_cols=inp_c
    print "false"
else:
    index = int(args.ind)
    inp_names1.append(inp_n[index])
    inp_names2=inp_n

print("Calculation Start")
print("index %s" % index)
sys.stdout.flush()
# Define Corr-Array
#corr_funcs=['msd','c1','c2'] # This is the name of the correlation functions.
corr_funcs=['c2']

# Read in input file
# Calculate important quantities
filesperblock=nfiles/float(nblocks)
t_val=stats.t.ppf(0.95,nblocks-1)/np.sqrt(nblocks)

# Creates the time array
time = np.genfromtxt('real_time.dat', usecols=0)
time = [0]+time
time = time[:ntimes]
for corr_name in corr_funcs:
    # Initialize Vectors
    cab=[ntimes]        # Unweighted Corr func
    fval1=[nfiles]      # Fluct Value 1
    fval2=[nfiles]      # Fluct Value 2
    fval3=[nfiles]      # Fluct Value 3
    fval4=[nfiles]      # Fluct Value 4
    dfval1=[nfiles]     # Fluct 1
    dfval2=[nfiles]     # Fluct 2
    dfval3=[nfiles]     # Fluct 3
    dfval4=[nfiles]     # Fluct 4
    d2fval=[nfiles]     # Second order weighted correlation function
    d3fval=[nfiles]     # Third order weighted correlation function
    d4fval=[nfiles]     # Fourth order weighted correlation function
    
    # First Derivatives
    dfvalcab=[ntimes]   # Weighted correlation function
    
    # Second Derivatives
    d2fvalcab=[ntimes]
    d2fvalavcab=[ntimes]

    # Third Derivative
    d3fvalcab=[ntimes]
    
    # Fourth Derivative
    d4fvalcab=[ntimes]
    
    # Block Initializations
    eafval_bl=[[] for x in range(0,nblocks)]
    cab_bl=[[] for x in range(0,nblocks)]
    dcabfval_bl=[[] for x in range(0,nblocks)]
    d2cabfval_bl=[[] for x in range(0,nblocks)]
    d2cabavfval_bl=[[] for x in range(0,nblocks)]
    d3cabfval_bl=[[] for x in range(0,nblocks)]
    d4cabfval_bl=[[] for x in range(0,nblocks)]

    for i in range(0,nfiles):
        dfval1.append(0)
        dfval2.append(0)
        dfval3.append(0)
        dfval4.append(0)
        d2fval.append(0)
        d3fval.append(0)
        d4fval.append(0)

    for i in range(0,ntimes-1):
        cab.append(0)
        dfvalcab.append(0)
        d2fvalcab.append(0)
        d2fvalavcab.append(0)
        d3fvalcab.append(0)
        d4fvalcab.append(0)

    # Read in File Names
    fnames=np.genfromtxt('file_names', dtype='string', unpack=True)
    fcab=["FILES/"+str(s)+"/"+str(corr_name)+"_"+str(s)+"_"+str(mol_name)+".dat" for s in fnames]
    item_count=0
    # Do the Block Average Calculation
    iindex=0
    jindex=0
    print("Read in fluctuations")
    sys.stdout.flush()
    for item1 in inp_names1:
        print item1
        sys.stdout.flush()
        fval1 = np.genfromtxt(str(item1)+'_init.out',unpack=True) # Reads in the values of the energy
        jindex=0
        for item2 in inp_names2:
            if jindex >= iindex:
                print item2
                item3=item1
                item4=item1

                fval2 = np.genfromtxt(str(item2)+'_init.out', unpack=True) # Reads in the values of the energy
                # HIGH ORDER DERIVATIVES
                # This only works for diagonal terms
                fval3 = fval1
                fval4 = fval1
                # END HIGH ORDER DERIVATIVES
                for block in range(0,nblocks):
                    start=int(block*filesperblock)
                    end=int((block+1)*filesperblock)
                    print start,end
                    fval1avg_bl = np.average(fval1[start:end])
                    fval2avg_bl = np.average(fval2[start:end])
                    fval3avg_bl = np.average(fval3[start:end])
                    fval4avg_bl = np.average(fval4[start:end])
                    for j in range(start,end):
                        dfval1[j]=fval1[j]-fval1avg_bl  # de1 = e1 - <e1>
                        dfval2[j]=fval2[j]-fval2avg_bl  # de2 = e2 - <e2>
                        dfval3[j]=fval3[j]-fval3avg_bl  # de3 = e3 - <e3>
                        dfval4[j]=fval4[j]-fval4avg_bl  # de4 = e4 - <e4>
                        d2fval[j]=dfval1[j]*dfval2[j]   # d2 = de1*de2 = (e-<e1>)(e-<e2>)
                        d3fval[j]=d2fval[j]*dfval3[j]   # d3 = de2
                        d4fval[j]=d3fval[j]*dfval4[j]   # d4 = d3*de3
                    d2fvalavg_bl = np.average(d2fval[start:end]) # <d2>
                    d3fvalavg_bl = np.average(d3fval[start:end])
                    d4fvalavg_bl = np.average(d4fval[start:end])
                    for j in range(0,ntimes):
                        cab[j]=0
                        dfvalcab[j]=0
                        d2fvalcab[j]=0
                        d2fvalavcab[j]=0
                        d3fvalcab=0
                        d4fvalcab=0
                    for j in range(start,end):
                        tcab=np.genfromtxt(fcab[j], usecols=1,unpack=True)
                        # Weight by Fluct
                        cab=cab+tcab # Unweighted correlation function
                        dfvalcab=dfvalcab+tcab*dfval1[j] # Singly weighted correlation function
                        d2fvalcab=d2fvalcab+tcab*d2fval[j] # Doubly weighted correlation function
                        d2fvalavcab=d2fvalavcab+tcab*d2fvalavg_bl # Weighted correlation function by avg
                        d3fvalcab=d3fvalcab+tcab*Third_Order_Weight(j,d3fval, d3fvalavg_bl)
                        d4fvalcab=d4fvalcab+tcab*Fourth_Order_Weight(j,d4fval, d4fvalavg_bl)

                    # Add Contributions from other correlation functions
                    d3fvalcab=d3fvalcab+3*dfvalcab*d2fvalavg_bl
                    d4fvalcab=d4fvalcab-6*np.subtract(d2fvalcab,d2fvalavcab)*d2fvalavg_bl + 4*dfvalcab*d3fvalavg_bl
                                                            
                    # Normalize by the number of files in the block
                    cab[:] = [x / float(filesperblock) for x in cab]  
                    dfvalcab[:] = [x / float(filesperblock) for x in dfvalcab]
                    d2fvalcab[:] = [x / float(filesperblock) for x in d2fvalcab]
                    d2fvalavcab[:] = [x / float(filesperblock) for x in d2fvalavcab]
                    d3fvalcab[:] = [x / float(filesperblock) for x in d3fvalcab]
                    d4fvalcab[:] = [x / float(filesperblock) for x in d4fvalcab]


                    # Calculate the Ratio of Weighted to Unweighted Corr Functions
                    eafval=dfvalcab[1:]/cab[1:]
                    # Save to list
                    eafval=np.insert(eafval,0,0)

                    eafval_bl[block] = list(eafval)
                    # Save correlation info to sep list
                    cab_bl[block]=list(cab) # block values for correlation function
                    dcabfval_bl[block]=list(dfvalcab)  # block values for weighted correlation function
                    d2cabfval_bl[block]=list(d2fvalcab) # block values for doubly weighted correlation function
                    d2cabavfval_bl[block]=list(d2fvalavcab) # block values for doubly averaged correlation function
                    d3cabfval_bl[block]=list(d3fvalcab) # block values for triply weighted correlation function
                    d4cabfval_bl[block]=list(d4fvalcab) # block values for quadruply weighted correlation function
                    # Print Single Correlation Values
                    np.savetxt("bl_"+str(block)+"_"+item1+"_"+mol_name+"_"+corr_name+".dat", np.c_[time, cab, dfvalcab,eafval], fmt='%s')
                    # Print Second Derivatives
                    np.savetxt("bl_"+str(block)+"_"+item1+"_"+item2+"_"+mol_name+"_"+corr_name+".dat", np.c_[time,d2fvalcab, d2fvalavcab, d2fvalcab-d2fvalavcab])
                    # Print Third Derivatives
                    np.savetxt("bl_"+str(block)+"_"+item1+"_"+item2+"_"+item3+"_"+mol_name+'_'+corr_name+".dat", np.c_[time, d3fvalcab])
                    # Print Fourth Derivatives
                    np.savetxt("bl_"+str(block)+"_"+item1+"_"+item2+"_"+item3+"_"+item4+"_"+mol_name+'_'+corr_name+".dat", np.c_[time, d4fvalcab])
                    # Zero the items
                    for i in range(start,end):
                        dfval1[i]=0
                        dfval2[i]=0
                        dfval3[i]=0
                        dfval4[i]=0
                        d2fval[i]=0
                        d3fval[i]=0
                        d4fval[i]=0
                    for i in range(0,ntimes):
                        cab[i]=0
                        dfvalcab[i]=0
                        d2fvalcab[i]=0
                        d2fvalavcab[i]=0
                        d3fvalcab[i]=0
                        d4fvalcab[i]=0

                    for i in range(0,ntimes):
                        eafval[i]=0

            err_cab = np.array(cab_bl).std(0)
            err_eafvalcab = np.array(eafval_bl).std(0)
            err_dfvalcab = np.array(dcabfval_bl).std(0)
            err_d2fvalcab = np.array(d2cabfval_bl).std(0)
            err_d2fvalavcab = np.array(d2cabavfval_bl).std(0)
            err_d3fvalcab = np.array(d3cabfval_bl).std(0)
            err_d4fvalcab = np.array(d4cabfval_bl).std(0)

            err_cab = [x * t_val for x in err_cab]
            err_eafvalcab = [x * t_val for x in err_eafvalcab]
            err_dfvalcab = [x * t_val for x in err_dfvalcab]
            err_d2fvalcab = [x * t_val for x in err_d2fvalcab]
            err_d2fvalavcab = [x * t_val for x in err_d2fvalavcab]
            err_d3fvalcab = [x * t_val for x in err_d3fvalcab]
            err_d4fvalcab = [x * t_val for x in err_d4fvalcab]

            print len(err_cab)
            print len(err_eafvalcab)
            print len(err_dfvalcab)
            np.savetxt('err_'+item1+'_'+str(mol_name)+'_'+corr_name+".dat", np.c_[err_cab, err_eafvalcab, err_dfvalcab], fmt='%s')
            np.savetxt('err_'+item1+'_'+item2+'_'+str(mol_name)+'_'+corr_name+'.dat', np.c_[err_d2fvalcab, err_d2fvalavcab], fmt='%s')
            np.savetxt('err_'+item1+'_'+item2+'_'+item3+'_'+item4+'_'+mol_name+'_'+corr_name+'.dat', np.c_[err_d3fvalcab, err_d4fvalcab], fmt='%s')
            jindex+=1
        item_count+=1
        iindex+=1
    # This Ends the Block Calculation
    for i in range(0,nfiles):
        fval1[i]=0
        fval2[i]=0
        fval3[i]=0
        fval4[i]=0
        dfval1[i]=0
        dfval2[i]=0
        dfval3[i]=0
        dfval4[i]=0
        d2fval[i]=0
        d3fval[i]=0
        d4fval[i]=0

    # Do the total calculation
    item_count = 0
    iindex=0 #item1 counter
    jindex=0 #item2 counter
    print("Total Calculation")
    for item1 in inp_names1:
        print("item1: %s" % item1)
        sys.stdout.flush()
        jindex=0
        fval1 = np.genfromtxt(str(item1)+"_init.out", unpack=True) # Reads in Values for the Energy
        for item2 in inp_names2:
            if jindex >= iindex:
                print("item2: %s" % item2)
                item3 = item1
                item4 = item1
                fval2 = np.genfromtxt(str(item2)+'_init.out', unpack=True) # Reads int he values of the energy
                # HIGH ORDER DERIVATIVES - DIAGONAL TERMS ONLY 
                # This would need to be expanded into a double loop to calc the full combinations.
                fval3 = fval1
                fval4 = fval1
                # END HIGH ORDER DERIVATIVES 
                fval1avg = np.average(fval1)
                fval2avg = np.average(fval2)
                fval3avg = np.average(fval3)
                fval4avg = np.average(fval4)
                print fval1avg, fval2avg
                for j in range(0, nfiles):
                    dfval1[j]=fval1[j]-fval1avg # de1 = e1 - <e1>
                    dfval2[j]=fval2[j]-fval2avg # de2 = e2 - <e2>
                    dfval3[j]=fval3[j]-fval3avg # de3 = e3 - <e3>
                    dfval4[j]=fval4[j]-fval4avg # de4 = e4 - <e4>
                    d2fval[j]=dfval1[j]*dfval2[j] # d2e = de1*de2 = (e-<e1>)*(e-<e2>)
                    d3fval[j]=d2fval[j]*dfval3[j] # d3e = d2e*de3
                    d4fval[j]=d3fval[j]*dfval4[j] # d4e = d3e*de4
                d1fvalavg = np.average(dfval1[:])
                d2fvalavg = np.average(d2fval[:]) # <d2e>
                d3fvalavg = np.average(d3fval[:])
                d4fvalavg = np.average(d4fval[:])
                print("Average Values:\n")
                print("<d1>: %s" % d1fvalavg)
                print("<d2>: %s" % d2fvalavg)
                print("<d3>: %s" % d3fvalavg)
                print("<d4>: %s" % d4fvalavg)
                # Zeroes all the correlation function arrays
                for j in range(0,ntimes):
                    cab[j]=0
                    dfvalcab[j]=0
                    d2fvalcab[j]=0
                    d2fvalavcab[j]=0
                    d3fvalcab[j]=0
                    d4fvalcab[j]=0
                for j in range(0,nfiles):
                    tcab=np.genfromtxt(fcab[j],usecols=1, unpack=True)
                    # Weight by flucts
                    cab=cab+tcab
                    dfvalcab=dfvalcab+tcab*dfval1[j] # Weighted Correlation Function
                    d2fvalcab=d2fvalcab+tcab*d2fval[j] # Doubly Weighted Correlation Function
                    d2fvalavcab=d2fvalavcab+tcab*d2fvalavg # Doubly weighted average correlation function
                    d3fvalcab=d3fvalcab+tcab*Third_Order_Weight(j,d3fval, d3fvalavg) # Triply weighted correlation function
                    d4fvalcab=d4fvalcab+tcab*Fourth_Order_Weight(j,d4fval, d4fvalavg) # Quadruply weighted correlation function
                # Calculated Lower Order Contributions
                d3fvalcab=d3fvalcab+3*dfvalcab*d2fvalavg
                d4fvalcab=d4fvalcab - 6*np.subtract(d2fvalcab,d2fvalavcab)*d2fvalavg - 4.0*dfvalcab*d3fvalavg
                # Normalize
                cab[:] = [x / float(nfiles) for x in cab]
                dfvalcab[:] = [x / float(nfiles) for x in dfvalcab]
                d2fvalcab[:] = [x / float(nfiles) for x in d2fvalcab]
                d2fvalavcab[:] = [x / float(nfiles) for x in d2fvalavcab]
                d3fvalcab[:] = [x / float(nfiles) for x in d3fvalcab]
                d4fvalcab[:] = [x / float(nfiles) for x in d4fvalcab]
                
                # Calculate the Ratio of Weighted to Unweighted Corr Functions
                eafval=dfvalcab[1:]/cab[1:]

                # Read In Uncertainty
                err_cab, err_eafvalcab, err_dfvalcab = np.genfromtxt('err_'+item1+'_'+str(mol_name)+'_'+corr_name+'.dat', usecols=(0,1,2), unpack=True)
                err_d2fvalcab, err_d2fvalavcab = np.genfromtxt('err_'+item1+'_'+item2+'_'+str(mol_name)+'_'+corr_name+'.dat',usecols=(0,1), unpack=True)
                err_d3fvalcab, err_d4fvalcab = np.genfromtxt('err_'+item1+'_'+item2+'_'+item3+'_'+item4+'_'+mol_name+'_'+corr_name+'.dat', usecols=(0,1), unpack=True)
                # Print Output
                eafval = np.insert(eafval,0,0)
                np.savetxt(item1+"_"+mol_name+"_"+corr_name+".dat", np.c_[time, cab, err_cab, dfvalcab,err_dfvalcab, eafval, err_eafvalcab], fmt='%s')
                np.savetxt(item1+"_"+item2+"_"+mol_name+"_"+corr_name+".dat", np.c_[time, d2fvalcab, err_d2fvalcab, d2fvalavcab,err_d2fvalavcab, np.subtract(d2fvalcab,d2fvalavcab), np.subtract(err_d2fvalcab,err_d2fvalavcab)], fmt='%s')
                np.savetxt(item1+"_"+item2+"_"+item3+"_"+mol_name+"_"+corr_name+".dat", np.c_[time, d3fvalcab, err_d3fvalcab], fmt='%s')
                np.savetxt(item1+"_"+item2+"_"+item3+"_"+item4+"_"+mol_name+"_"+corr_name+".dat", np.c_[time, d4fvalcab, err_d4fvalcab], fmt='%s')
                # Zero Items
                for i in range(0,nfiles):
                    dfval1[i]=0
                    dfval2[i]=0
                    dfval3[i]=0
                    dfval4[i]=0
                    d2fval[i]=0
                    d3fval[i]=0
                    d4fval[i]=0
                    
                for i in range(0,ntimes):
                    cab[i]=0
                    dfvalcab[i]=0
                    d2fvalcab[i]=0
                    d2fvalavcab[i]=0
                    d3fvalcab[i]=0
                    d4fvalcab[i]=0

                for i in range(0,ntimes):
                    eafval[i]=0
            jindex+=1
        item_count += 1
        iindex+=1

# Program Complete
