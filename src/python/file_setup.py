#!/usr/bin/python

import sys
import os
import numpy as np
import math
from read_input import input

inputparam = input('input_file')

f = open('file_names','r')
ext = f.readlines()
g = open('mol_names', 'r')
molnames = []
molnames = g.readlines()

ext = [x.strip() for x in ext]
molnames = [x.strip() for x in molnames]

nfiles = len(ext)
nmols = len(molnames)
narrays=int(math.ceil(float(nfiles)/500.))

if inputparam.prog == "LAMMPS":
    narrays=int(math.ceil(float(nfiles)/500.))

    m = open ('sub_script', 'w')
    for i in range(0,int(narrays)):
        nstart=i*500+1
        nend=i*500+500
        m.write("cp job_array.sh job_array_"+str(i)+"\n")
        m.write("sed -i -e 's@AAA@"+str(int(nstart))+"@g' job_array_"+str(i)+"\n")
        m.write("sed -i -e 's@BBB@"+str(int(nend))+"@g' job_array_"+str(i)+"\n")
        m.write("msub job_array_"+str(i)+"\n")
        if i != narrays-1:
            if (nend % 5000) == 0:
                m.write("sleep 30m\n")
    m.write('touch .flag_nvecomplete')
elif inputparam.prog == "CP2K":
    narrays=int(math.ceil(float(nfiles)/500.))

    m = open ('sub_script', 'w')
    for i in range(0,int(narrays)):
        nstart=i*500+1
        nend=i*500+500
        m.write("cp job_array.sh job_array_"+str(i)+"\n")
        m.write("sed -i -e 's@AAA@"+str(int(nstart))+"@g' job_array_"+str(i)+"\n")
        m.write("sed -i -e 's@BBB@"+str(int(nend))+"@g' job_array_"+str(i)+"\n")
        m.write("msub job_array_"+str(i)+"\n")
        if i != narrays-1:
            if (nend % 5000) == 0:
                m.write("sleep 30m\n")
    m.write('touch .flag_nvecomplete')
else:  
    print("Input File Choice of %s is incorrect" % inputparam.prog)
           
