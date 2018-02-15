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
inp_names, inp_cols=np.genfromtxt(inputfile, usecols=(0,1), dtype=(str,int),unpack=True)
print inp_names
print inp_cols
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
    dfval1=[nfiles]     # Fluct 1
    dfval2=[nfiles]     # Fluct 2
    d2fval=[nfiles]     # Second order weighted correlation function
    
    # First Derivatives
    dfvalcab=[ntimes]   # Weighted correlation function
    
    # Second Derivatives
    d2fvalcab=[ntimes]
    d2fvalavcab=[ntimes]
    
    # Block Initializations
    eafval_bl=[[] for x in range(0,nblocks)]
    cab_bl=[[] for x in range(0,nblocks)]
    dcabfval_bl=[[] for x in range(0,nblocks)]
    d2cabfval_bl=[[] for x in range(0,nblocks)]
    d2cabavfval_bl=[[] for x in range(0,nblocks)]

    for i in range(0,nfiles):
        dfval1.append(0)
        dfval2.append(0)
        d2fval.append(0)

    for i in range(0,ntimes-1):
        cab.append(0)
        dfvalcab.append(0)
        d2fvalcab.append(0)
        d2fvalavcab.append(0)
    print len(cab)
    # Read in File Names
    fnames=np.genfromtxt('file_names', dtype='string', unpack=True)
    fcab=["FILES/"+str(s)+"/"+str(corr_name)+"_"+str(s)+"_"+str(mol_name)+".dat" for s in fnames]
    item_count=0
    # Do the Block Average Calculation
    iindex=0
    jindex=0
    for item1 in inp_names:
        print item1
        fval1 = np.genfromtxt(str(item1)+'_init.out',unpack=True) # Reads in the values of the energy
        jindex=0
        for item2 in inp_names:
            if jindex >= iindex:
                print item2
                fval2 = np.genfromtxt(str(item2)+'_init.out', unpack=True) # Reads in the values of the energy
                for block in range(0,nblocks):
                    start=int(block*filesperblock)
                    end=int((block+1)*filesperblock)
                    print start,end
                    fval1avg_bl = np.average(fval1[start:end])
                    fval2avg_bl = np.average(fval2[start:end])
                    for j in range(start,end):
                        dfval1[j]=fval1[j]-fval1avg_bl  # de1 = e1 - <e1>
                        dfval2[j]=fval2[j]-fval2avg_bl  # de2 = e2 - <e2>
                        d2fval[j]=dfval1[j]*dfval2[j]   # de2 = de1*de2 = (e-<e1>)(e-<e2>)
                    d2fvalavg_bl = np.average(d2fval[start:end]) # <de2>
                    for j in range(0,ntimes):
                        cab[j]=0
                        dfvalcab[j]=0
                        d2fvalcab[j]=0
                        d2fvalavcab[j]=0
                    for j in range(start,end):
                        tcab=np.genfromtxt(fcab[j], usecols=1,unpack=True)
                        # Weight by Fluct
                        cab=cab+tcab # Unweighted correlation function
                        dfvalcab=dfvalcab+tcab*dfval1[j] # Singly weighted correlation function
                        d2fvalcab=d2fvalcab+tcab*d2fval[j] # Doubly weighted correlation function
                        d2fvalavcab=d2fvalavcab+tcab*d2fvalavg_bl # Weighted correlation function by avg
                    # Normalize by the number of files in the block
                    cab[:] = [x / float(filesperblock) for x in cab]  
                    dfvalcab[:] = [x / float(filesperblock) for x in dfvalcab]
                    d2fvalcab[:] = [x / float(filesperblock) for x in d2fvalcab]
                    d2fvalavcab[:] = [x / float(filesperblock) for x in d2fvalavcab]
                    # Calculate the Ratio of Weighted to Unweighted Corr Functions
                    eafval=dfvalcab[1:]/cab[1:]
                    # Save to list
                    eafval_bl[block] = list(eafval)
                    # Save correlation info to sep list
                    cab_bl[block]=list(cab) # block values for correlation function
                    dcabfval_bl[block]=list(dfvalcab[1:])  # block values for weighted correlation function
                    d2cabfval_bl[block]=list(d2fvalcab[1:]) # block values for doubly weighted correlation function
                    d2cabavfval_bl[block]=list(d2fvalavcab[1:]) # block values for doubly averaged correlation function
                    # Print Single Correlation Values
                    np.savetxt("bl_"+str(block)+"_"+item1+"_"+mol_name+"_"+corr_name+".dat", np.c_[time[1:], cab[1:], dfvalcab[1:],eafval], fmt='%s')
                    # Print Second Derivatives
                    np.savetxt("bl_"+str(block)+"_"+item1+"_"+item2+"_"+mol_name+"_"+corr_name+".dat", np.c_[time[1:],d2fvalcab[1:], d2fvalavcab[1:], d2fvalcab[1:]-d2fvalavcab[1:]])
                    # Zero the items
                    for i in range(start,end):
                        dfval1[i]=0
                        dfval2[i]=0
                        d2fval[i]=0
                    for i in range(0,ntimes):
                        cab[i]=0
                        dfvalcab[i]=0
                        d2fvalcab[i]=0
                        d2fvalavcab[i]=0
                    for i in range(0,ntimes-1):
                        eafval[i]=0

            err_cab = np.array(cab_bl).std(0)
            err_eafvalcab = np.array(eafval_bl).std(0)
            err_dfvalcab = np.array(dcabfval_bl).std(0)
            err_d2fvalcab = np.array(d2cabfval_bl).std(0)
            err_d2fvalavcab = np.array(d2cabavfval_bl).std(0)

            err_cab = [x * t_val for x in err_cab]
            err_eafvalcab = [x * t_val for x in err_eafvalcab]
            err_dfvalcab = [x * t_val for x in err_dfvalcab]
            err_d2fvalcab = [x * t_val for x in err_d2fvalcab]
            err_d2fvalavcab = [x * t_val for x in err_d2fvalavcab]
            np.savetxt('err_'+item1+'_'+str(mol_name)+'_'+corr_name+".dat", np.c_[err_cab[1:], err_eafvalcab, err_dfvalcab], fmt='%s')
            np.savetxt('err_'+item1+'_'+item2+'_'+str(mol_name)+'_'+corr_name+'.dat', np.c_[err_d2fvalcab, err_d2fvalavcab], fmt='%s')
            jindex+=1
        item_count+=1
        iindex+=1
    # This Ends the Block Calculation
    for i in range(0,nfiles):
        fval1[i]=0
        fval2[i]=0
        dfval1[i]=0
        dfval2[i]=0
        d2fval[i]=0
    # Do the total calculation
    item_count = 0
    iindex=0 #item1 counter
    jindex=0 #item2 counter
    for item1 in inp_names:
        print item1
        jindex=0
        fval1 = np.genfromtxt(str(item1)+"_init.out", unpack=True) # Reads in Values for the Energy
        for item2 in inp_names:
            if jindex >= iindex:
                print item2
                fval2 = np.genfromtxt(str(item2)+'_init.out', unpack=True) # Reads int he values of the energy
                fval1avg = np.average(fval1)
                fval2avg = np.average(fval2)
                print fval1avg, fval2avg
                for j in range(0, nfiles):
                    dfval1[j]=fval1[j]-fval1avg # de1 = e1 - <e1>
                    dfval2[j]=fval2[j]-fval2avg # de2 = e2 - <e2>
                    d2fval[j]=dfval1[j]*dfval2[j] # d2e = de1*de2 = (e-<e1>)*(e-<e2>)
                d2fvalavg = np.average(d2fval[:]) # <d2e>
                for j in range(0,ntimes):
                    cab[j]=0
                    dfvalcab[j]=0
                    d2fvalcab[j]=0
                    d2fvalavcab[j]=0
                for j in range(0,nfiles):
                    tcab=np.genfromtxt(fcab[j],usecols=1, unpack=True)
                    # Weight by flucts
                    cab=cab+tcab
                    dfvalcab=dfvalcab+tcab*dfval1[j] # Weighted Correlation Function
                    d2fvalcab=d2fvalcab+tcab*d2fval[j] # Doubly Weighted Correlation Function
                    d2fvalavcab=d2fvalavcab+tcab*d2fvalavg # Doubly weighted average correlation function
                # Normalize
                cab[:] = [x / float(nfiles) for x in cab]
                dfvalcab[:] = [x / float(nfiles) for x in dfvalcab]
                d2fvalcab[:] = [x / float(nfiles) for x in d2fvalcab]
                d2fvalavcab[:] = [x / float(nfiles) for x in d2fvalavcab]
                
                # Calculate the Ratio of Weighted to Unweighted Corr Functions
                eafval=dfvalcab[1:]/cab[1:]

                # Read In Uncertainty
                err_cab, err_eafval, err_dfvalcab = np.genfromtxt('err_'+item1+'_'+str(mol_name)+'_'+corr_name+'.dat', usecols=(0,1,2), unpack=True)
                err_d2fvalcab, err_d2fvalavcab = np.genfromtxt('err_'+item1+'_'+item2+'_'+str(mol_name)+'_'+corr_name+'.dat',usecols=(0,1), unpack=True)
                # Print Output
                np.savetxt(item1+"_"+mol_name+"_"+corr_name+".dat", np.c_[time[1:], cab[1:], err_cab, dfvalcab[1:],err_dfvalcab, eafval, err_eafvalcab], fmt='%s')
                np.savetxt(item1+"_"+item2+"_"+mol_name+"_"+corr_name+".dat", np.c_[time[1:], d2fvalcab[1:], err_d2fvalcab, d2fvalavcab[1:],err_d2fvalavcab, np.subtract(d2fvalcab,d2fvalavcab)[1:], np.subtract(err_d2fvalcab,err_d2fvalavcab)], fmt='%s')
                # Zero Items
                for i in range(0,nfiles):
                    dfval1[i]=0
                    dfval2[i]=0
                    d2fval[i]=0
                for i in range(0,ntimes):
                    cab[i]=0
                    dfvalcab[i]=0
                    d2fvalcab[i]=0
                    d2fvalavcab[i]=0
                for i in range(0,ntimes-1):
                    eafval[i]=0
            jindex+=1
        item_count += 1
        iindex+=1

# Program Complete
