#!/usr/bin/env python
""" 
I realize this isn't a fortran program - but it was much simpler to write.
"""
import numpy as np

def min_dist(a, b):
    dr      = []
    shift   = []
    dsq     = 0.0
    for i in range(3):
        shift.append(L*round((a[i]-b[i])/L))
        dr.append(a[i]-b[i]-shift[i])
        dsq += dr[i]**2.
    drfinal = np.sqrt(dsq)
    return drfinal

def calc_ang(dOO,dHO):
    dHOO = np.arccos(np.around((0.9572**2. + dOO**2. - dHO**2.)/(2.*0.9572*dOO),4))
    return dHOO
    
def read_frames(filename):
    rO = []
    r1 = []
    r2 = []
    
    with open(filename) as f:
        count = 0
        nmols = 0
        index = 0
        for line in f:
            if len(line.split()) == 1:
                count = 0
                nmols = 0
            if len(line.split()) == 4:
                # Check if oxygen
                if count % 3 == 0:
                    # Read oxygen atoms and wrap boundary conditions
                    rO.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
                elif count % 3 == 1:
                    # Read in hydrogens and wrap boundary conditions
                    r1x = float(line.split()[1]) - L*round((rO[index][0]-float(line.split()[1]))/L)
                    r1y = float(line.split()[2]) - L*round((rO[index][1]-float(line.split()[2]))/L)
                    r1z = float(line.split()[3]) - L*round((rO[index][1]-float(line.split()[3]))/L)
                    r1.append([r1x, r1y, r1z])
                elif count % 3 == 2:
                    # Read in hydrogens and wrap boundary conditions
                    r2x = float(line.split()[1]) - L*round((rO[index][0]-float(line.split()[1]))/L)
                    r2y = float(line.split()[2]) - L*round((rO[index][1]-float(line.split()[2]))/L)
                    r2z = float(line.split()[3]) - L*round((rO[index][2]-float(line.split()[3]))/L)
                    r2.append([r2x, r2y, r2z])
                    index += 1
                    nmols += 1
                else:
                    print("Incorrect count during read")
                count += 1
    return nmols, rO, r1, r2
# Read Corr Calc Input File
inpfile='corr_calc.in'
inp = open(inpfile, 'r')
lines = []
for line in inp:
    lines.append(line)
nfile  = lines[1]
ntimes = int(lines[3].split()[0])
dt     = float(lines[3].split()[1])
volume = float(lines[5])
mol_name = lines[7]
inp.close()

L = volume**(1/3.)
rOO_max = 3.1
rHO_max = 2.0
ang_max = 20.

filename = 'traj.xyz'

# Read the frames
nmols, rO, r1, r2 = read_frames(filename)

# Determine initial HBONDS

OHs =[]
for mol1 in range(nmols):
    for mol2 in range(nmols):
        if mol1 != mol2:
            dOO = min_dist(rO[mol1],rO[mol2])
            if dOO < rOO_max:
                dHO1 = min_dist(r1[mol1], rO[mol2])
                dHO2 = min_dist(r2[mol1], rO[mol2])
                if dHO1 < rHO_max:
                    dHOO1 = np.degrees(calc_ang(dOO, dHO1))
                    if dHOO1 < ang_max:
                        OHs.append([mol1,mol2,1,])
                if dHO2 < rHO_max:
                    dHOO2 = np.degrees(calc_ang(dOO, dHO2))
                    if dHOO2 < ang_max:
                        OHs.append([mol1,mol2,2])

print("There are %s OHs" % len(OHs))

crp = np.zeros(ntimes)
steps = [0]
for n in range(1,ntimes):
    crp[n] = crp[n-1]
    steps.append(n)
    timeindex = n*nmols
    for OH in range(len(OHs)):
        mol1 = OHs[OH][0]
        for mol2 in range(nmols):
            if OHs[OH][0] != mol2 and OHs[OH][2] != 0:
                dOO = min_dist(rO[mol1+timeindex],rO[mol2+timeindex])
                assert dOO != 0, "dOO = 0... NO!"
                if dOO < rOO_max:
                    dHO1 = min_dist(r1[mol1+timeindex], rO[mol2+timeindex])
                    dHO2 = min_dist(r2[mol1+timeindex], rO[mol2+timeindex])
                    if dHO1 < rHO_max and OHs[OH][2] == 1:
                        dHOO1 = np.degrees(calc_ang(dOO, dHO1))
                        if dHOO1 < ang_max:
                            if OHs[OH][1] == mol2:
                                crp[n]+=0
                            else:
                                crp[n]+=1
                                OHs[OH][2]=0
                    if dHO2 < rHO_max and OHs[OH][2] == 2:
                        dHOO2 = np.degrees(calc_ang(dOO, dHO2))
                        if dHOO2 < ang_max:
                            if OHs[OH][1] == mol2:
                                crp[n]+=0
                            else:
                                crp[n]+=1
                                OHs[OH][2]=0

np.savetxt('test.out', np.c_[steps, crp/float(len(OHs))], fmt="%2.5f")



                
