import numpy as np

it=0
with open ("log.lammps","r") as log:
    for line in log:
        it+=1
        if 'Pxx Pyy Pzz' in line:
            start=it
            
        if 'Loop time' in line:
            end=it

P = np.genfromtxt('log.lammps', skip_header=start, skip_footer=it-end, usecols=(13,14,15,16,17,18))

np.savetxt('pressures_out.log', P, fmt='%f8')


print start, end