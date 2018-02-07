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
parser.add_argument('-ntimes', help="Number of times")
parser.add_argument('-mol', help ="Molecule Name")
args = parser.parse_args()

inputfile = str(args.inp)
nfiles = int(args.files)
nblocks = int(args.blocks)
ntimes = int(args.ntimes)
mol_name = str(args.mol)

# Define Corr-Array
corr_funcs=['msd','c1','c2']

# Read in input file
inp_names, inp_cols=np.genfromtxt(inputfile, usecols=(0,1), dtype=None,unpack=True)

# Calculate important quantities
filesperblock=nfiles/float(nblocks)
t_val=stats.t.ppf(0.95,nblocks-1)/np.sqrt(nblocks)

time = []
for i in range(0, ntimes):
    time.append(i*0.05)

for corr_name in corr_funcs:
    # Initialize Vectors
    cab=[ntimes]        # Unweighted Corr func
    fval1=[nfiles]      # Fluct Value 1
    fval2=[nfiles]      # Fluct Value 2
    dfval1=[nfiles]     # Weighted Fluct Correlation Function 1
    dfval2=[nfiles]     # Weighted Fluct Correlation Function 2
    
    # First Derivatives
    dfvalcab=[ntimes]
    
    # Second Derivatives
    d2fvalcab=[ntimes]
    
    # Block Initializations
    eafval_bl[[] for x in range(0,nblocks)]
    cab_bl[[] for x in range(0,nblocks)]
    dcabfval_bl=[[] for x in range(0,nblocks)]
    
    for i in range(0,nfiles):
        dfval1.append(0)
        dfval2.append(0)

    # Read in File Names
    fnames=[np.genfromtxt('file_names', dtype='string', unpack=True)]
    fcab=["FILES/"+s+"/"+corr_name+"_"+s+"_"+mol_name+".dat" for s in fnames]
    item_count=0
    # Do the Block Average Calculation
    for item1 in inp_names:
        print item1
        fval1 = np.genfromtxt(item1+'_init.out',unpack=True) # Reads in the values of the energy
        for item2 in inp_names:
            if inp_names.index[item2] >= item_count:
                print item2
                fval2 = np.genfromtxt(item2+'_init.out', unpack=True) # Reads in the values of the energy
                for block in range(0,nblocks):
                    start=int(block*filesperblock)
                    end=int((block+1)*filesperblock)
                    print start,end
                    fval1avg_bl = np.average(fval1[start:end])
                    fval2avg_bl = np.average(fval2[start:end])
                    for j in range(start,end):
                        dfval1[j]=fval1[j]-fval1avg_bl
                        dfval2[j]=fval2[j]-fval2avg_bl
                        d2fval[j]=dfval1[j]*dfval2[j]
                    d2fvalavg_bl = np.average(d2fval)
                    for j in range(0,ntimes):
                        cab[j]=0
                        dfvalcab[j]=0
                        d2fvalcab[j]=0
                        d2fvalavcab[j]=0
                    for j in range(start,end):
                        tcab=np.genfromtxt(fcab[j], usecols=1,unpack=True)
                        # Weight by Fluct
                        cab=cab+tcab
                        dfvalcab=dfvalcab+tcab*dfval[j]
                        d2fvalcab=d2fvalcab+tcab*d2fval[j]
                        d2fvalavcab=d2fvalavcab+tcab*d2fvalavg_bl
                    # Normalize
                    cab[:] = [x / float(filesperblock) for x in cab]
                    dfvalcab[:] = [x / float(filesperblock) for x in dfvalcab]
                    d2fvalcab[:] = [x / float(filesperblock) for x in d2fvalcab]
                    d2fvalavcab[:] = [x / float(filesperblock) for x in d2fvalavcab]
                    # Calculate the Ratio of Weighted to Unweighted Corr Functions
                    eafval=dfval[1:]/cab[1:]
                    # Save to list
                    eafval_bl[block] = list(eafval)
                    # Save correlation info to sep list
                    cab_bl[block]=list(cab)
                    dcabfval_bl[block]=list(dfvalcab[1:])
                    d2cabfval_bl[block]=list(d2fvalcab[1:])
                    d2cabavfval_bl[block]=list(d2fvalavcab[1:])
                    # Print Single Correlation Values
                    np.savetxt("bl_"+str(block)+"_"+item+"_"+mol_name+"_"+corr_name+".dat", np.c_[time[1:], cab[1:], dfvalcab[1:],eafval], fmt='%s')
                    # Print Second Derivatives
                    np.savetxt("bl_"+str(block)+"_"+item1+"_"+item2+"_"+mol_name+"_"+corr_name+".dat", np.c_[time[1:],d2fvalcab[1:], d2fvalavcab[1:], d2fvalcab[1:]-d2fvalavcab[1:])
                    # Zero the items
                    for i in range(start,end):
                        dfval1[i]=0
                        dfval2[i]=0
                    for i in range(0,ntimes):
                        cab[i]=0
                        dfvalcab[i]=0
                        d2fvalcab[i]=0
                        d2fvalavcab[i]=0
                    for i in range(0,ntimes-1):
                        eafval[i]=0
        item_count+=1
    # This Ends the Block Calculation
    # Do the total calculation
    item_count = 0
    for item1 in inp_names:
        print item1
        fval1 = np.genfromtxt(item1+"_init.out", unpack=True) # Reads in Values for the Energy
        for item2 in inp_names:
            if inp_names.index[item2] >= item_count:
                print item 2
                fval2 = np.genfromtxt(item2+'_init.out', unpack=True) # Reads int he values of the energy
                fval1avg = np.average(fval1)
                fval2avg = np.average(fval2)
                for i in range(0, nfiles):
                    dfval1[j]=fval1[j]-fval1avg
                    dfval2[j]=fval2[j]-fval2avg
                    d2fval[j]=dfval1[j]*dfval2[j]
                d2fvalavg = np.average(d2fval)
                for j in range(0,ntimes):
                    cab[j]=0
                    dfvalcab[j]=0
                    d2fvalcab[j]=0
                    d2fvalavcab[j]=0
                for j in range(0,nfiles):
                    tcab=np.genfromtxt(fcab[j],usecols=1, unpack=True)
                    # Weight by flucts
                    cab=cab+tcab
                    dfvalcab=dfvalcab+tcab*dfval[j]
                    d2fvalcab=d2fvalcab+tcab*d2fval[j]
                    d2fvalavcab=d2fvalavcab+tcab*d2fvalavg
                # Normalize
                cab[:] = [x / float(nfiles) for x in cab]
                dfvalcab[:] = [x / float(nfiles) for x in dfvalcab]
                d2fvalcab[:] = [x / float(nfiles) for x in d2fvalcab]
                d2fvalavcab[:] = [x / float(nfiles) for x in d2fvalavcab]
                # Calculate the Ratio of Weighted to Unweighted Corr Functions
                eafval=dfval[1:]/cab[1:]

                # Calculate Uncertainty
                err_cab = np.array(cab_bl).std(0)
                err_fvalcab = np.array(eafval_bl).std(0)
                err_dfvalcab = np.array(dfvalcab_bl).std(0)
                err_d2fvalcab = np.array(dfvalcab_bl).std(0)
                err_d2fvalavcab = np.array(dfvalavcab_bl).std(0)
                
                err_cab = [x * t_val for x in err_cab]
                err_fvalcab = [x * t_val for x in err_fvalcab]
                err_dfvalcab = [x * t_val for x in err_dfvalcab]
                err_d2fvalcab = [x * t_val for x in err_d2fvalcab]
                err_d2fvalavcab = [x * t_val for x in err_d2fvalavcab]

                # Print Output
                np.savetxt(item1+"_"+mol_name+"_"+corr_name+".dat", np.c_[time[1:], cab[1:], err_cab[1:], dfvalcab[1:],err_dfvalcab[1:], eafval,err_fvalcab], fmt='%s')
                np.savetxt(item1+"_"+item_2+"_"+mol_name+"_"+corr_name+".dat", np.c_[time[1:], d2fvalcab[1:], err_d2fvalcab[1:], d2fvalavcab[1:],err_d2fvalavcab[1:], d2fvalcab[1:]-d2fvalavcab[1:], err_d2fvalcab[1:]-err_d2fvalavcab[1:]], fmt='%s')
                # Zero Items
                for i in range(0,nfiles):
                    dfval1[i]=0
                    dfval2[i]=0
                for i in range(0,ntimes):
                    cab[i]=0
                    dfvalcab[i]=0
                    d2fvalcab[i]=0
                    d2fvalavcab[i]=0
                for i in range(0,ntimes-1):
                    eafval[i]=0
        item_count += 1

# Program Complete
