#!/usr/bin/python

import sys
import os
import numpy as np
import math

f = open('file_names','r')
ext = f.readlines()
g = open('mol_names', 'r')
molnames = []
molnames = g.readlines()

ext = [x.strip() for x in ext]
molnames = [x.strip() for x in molnames]

nfiles = len(ext)
nmols = len(molnames)

h = open('setup_files','w')

for i in range (0,nfiles):
    h.write('mkdir FILES/'+ext[i]+'\n')
    h.write('cp in.nve FILES/'+ext[i]+'\n')
    h.write('cp nve.sh FILES/'+ext[i]+'\n')
    h.write('cp set_msd_calcs.py FILES/'+ext[i]+'\n')
    h.write('cp msd_rot_calc FILES/'+ext[i]+'\n')
    h.write('cp grab_press.py FILES/'+ext[i]+'\n')
    h.write('cp visc_calc FILES/'+ext[i]+'\n')
    h.write('cp RESTART/restart.'+ext[i]+' FILES/'+ext[i]+'\n')
    h.write('cd FILES/'+ext[i]+'\n')
    h.write("sed -i -e 's@direct_calc_nve@nve_"+ext[i]+"@g' nve.sh\n")
    h.write("sed -i -e 's@restart.file@../../RESTART/restart."+ext[i]+"@g' in.nve\n")
    for j in range(0, nmols):
        h.write("sed -i -e 's@traj.file"+str(j)+"@traj_"+ext[i]+"_"+molnames[j]+".xyz@g' in.nve\n")
    h.write('cd ../../\n')

narrays=int(math.ceil(float(nfiles)/500.))

m = open ('sub_script', 'w')
for i in range(0,int(narrays)):
    nstart=i*500+1
    nend=i*500+500
    m.write("cp job_array.sh job_array_"+str(i)+"\n")
    m.write("sed -i -e 's@AAA@"+str(int(nstart))+"@g' job_array_"+str(i)+"\n")
    m.write("sed -i -e 's@BBB@"+str(int(nend))+"@g' job_array_"+str(i)+"\n")
    m.write("msub job_array_"+str(i)+"\n")
    if (nend % 5000) == 0:
        m.write("sleep 2h\n")

    
    
