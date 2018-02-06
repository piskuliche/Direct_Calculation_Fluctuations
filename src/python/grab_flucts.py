#grab_flucts.py is the program that parses each log file and grabs each of the energies/volumes/whatever thermodynamic data you want to weight the correlation functions
#required input: log.lammps (for each directory)
#produced ouptut: *_init.out 

import numpy as np

filenames='file_names'

input_file='flucts.inp'
avg=[]
fnames=[]



# Read Input File
ivalname, icolnum = np.genfromtxt(input_file, usecols=(0,1), unpack=True)

# Loop over types of values, e, v, ke, etc
for value in range(0, len(ivalname))
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
            ival = np.genfromtxt(filename, skip_header=startskip, skip_footer=endskip, usecols=(ivalnum-1), unpack=True)
            ival_list = np.append(ival_list, ival[0])

    fileout=ivalname+"_init.out"
    np.savetxt(fileout, ival_list, fmt=['%.4f'])

        
