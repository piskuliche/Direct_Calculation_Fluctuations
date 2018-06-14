#grab_flucts.py is the program that parses each log file and grabs each of the energies/volumes/whatever thermodynamic data you want to weight the correlation fucntions. 
#   required input: either log.lammps or flucts.ener
#   produced output: *_init.out

import numpy as np
import sys
from read_input import input

value = int(sys.argv[1])

inputparam = input('input_file')

filenames='file_names'

ivalname, icolnum = np.genfromtxt(input_file, usecols(0,1), dtype=(str,int), unpack=True)
ival_list = np.array([])

if inputparam.prog == 'LAMMPS':
    # Looks for log.lammps
    with open(filenames) as f:
        for l in f:
            filename='FILES/'+l.rstrip()+'/log.lammps'
            lookup='Step'
            lookup2='Loop time'
            totlines = 0
            with open(filename) as myFile:
                for num, line in enumerate(myFile,1):
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
else if inputparam.prog == 'CP2K':
    # Looks for CP2K Output File
    with open(filenames) as f:
        for l in f:
            filename='FILES/'+l.rstrip()+'flucts.ener'
            with open(filename) as myFile:
                myFile.readline()
                ival=(myFile.readline().strip())[int(icolnum[value])-1], unpack=True)
                ival_list = np.append(ival_list, ival)
    fileout=str(ivalname[value])+"_init.out"
    np.savetxt(fileout, ival_list, fmt=['%.4f'])

            
else:
    print("%s is not a valid program" % inputparam.prog)

