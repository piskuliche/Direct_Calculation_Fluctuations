#grab_flucts.py is the program that parses each log file and grabs each of the energies/volumes/whatever thermodynamic data you want to weight the correlation functions
#required input: log.lammps (for each directory)
#produced ouptut: *_init.out 

import numpy as np
import argparse
from argparse import RawTextHelpFormatter

# Read in Arguments
parser = argparse.ArgumentParser(description='''Calculates the fluctuations of the NPT Ensemble''', formatter_class=RawTextHelpFormatter)
parser.add_argument('-inp', help="Input File")
parser.add_argument('-ind', help = "If you want a single element - type its position in the input file (counting from 0). If you want the whole thing - type FALSE")
args = parser.parse_args()

filenames='file_names'

input_file=str(args.inp)
value=int(args.ind)
avg=[]
fnames=[]



# Read Input File
ivalname, icolnum = np.genfromtxt(input_file, usecols=(0,1), dtype=(str,int), unpack=True)

print ivalname[0]
print icolnum[0]

# Loop over types of values, e, v, ke, etc
ival_list = np.array([])
with open(filenames) as f:
    for l in f:
        filename='FILES/'+l.rstrip()+'/log.lammps'
        print filename
        lookup = 'Step'
        lookup2 = 'Loop time'
        totlines = 0
        with open(filename) as myFile:
            for num, line in enumerate(myFile, 1):
                if lookup in line:
                    startskip=num
                if lookup2 in line:
                    endskip=num
                totlines=num
        endskip=totlines-endskip-4
        ival = np.genfromtxt(filename, skip_header=startskip, skip_footer=endskip, usecols=(int(icolnum[value])-1), unpack=True)
        ival_list = np.append(ival_list, ival[0])

fileout=str(ivalname[value])+"_init.out"
np.savetxt(fileout, ival_list, fmt=['%.4f'])

        
