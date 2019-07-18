#!/usr/bin/env python
"""
This is a python code to calculate specific interactions and tabulate their interactions
    General Outline:
        1) read input file for pair coeffs
        2) read data file for types, molecules, charges and box length
        3) read traj file for xyz coordinates
        4) calculate interaction matrix
        5) caculate lennard jones potential
        6) calculate short range coulomb
"""
def read_paircoeffs(file):
    tmptype = []
    tmpeps =[]
    tmpsig = []
    with open(file,'r') as f:
        lines = f.readlines()
        for line in lines:
            if "pair_coeff" in line:
                tmptype.append(int(line.split()[1]))
                tmpeps.append(float(line.split()[3]))
                tmpsig.append(float(line.split()[4]))
    eps = []
    sig = []
    for ty in tmptype:
        eps.append(0)
        sig.append(0)
    for i in range(len(tmptype)):
        eps[tmptype[i]-1] = tmpeps[i]
        sig[tmptype[i]-1] = tmpsig[i]
    return eps, sig

def read_data(file):
    count = 0
    l_row = 0
    a_row = 0
    with open(file,'r') as f:
        lines = f.readlines()
        for line in lines:
            count += 1
            if "Masses" in line:
                l_stop = count - 1
                l_start = l_stop - 4
            if "Atoms" in line:
                a_start = count + 1 
            if "Bonds" in line:
                a_stop = count - 1
        l_row = int(l_stop - l_start - 1)
        a_row = int(a_stop - a_start - 1)
    llow,lhigh = np.genfromtxt(file,skip_header=l_start,max_rows=l_row,usecols=(0,1),unpack=True)
    tmpatom,tmpmol,tmptype,tmpq = np.genfromtxt(file,skip_header=a_start,max_rows=a_row,usecols=(0,1,2,3),unpack=True)
    mol = []
    typ = []
    Q = []
    L = []
    mol_count = int(tmpmol[-1])
    for c in range(mol_count):
        mol.append([])
    for a in range(len(tmpmol)):
        m = int(tmpmol[a] - 1)
        mol[m].append(int(tmpatom[a] - 1))
    for t in range(len(tmptype)):
        typ.append(int(tmptype[t]))
    for q in range(len(tmpq)):
        Q.append(tmpq[q])
    for l in range(len(llow)):
        L.append(lhigh[l]-llow[l])
    return mol, typ, Q, L

def read_frame(file):
    count = 0
    rx = []
    ry = []
    rz = []
    with open(file,'r') as f:
        lines = f.readlines()
        for line in lines:
            count += 1
            if "Atoms" in line:
                start = count 
        stop = count
        row = int(stop - start)
    rx, ry, rz = np.genfromtxt(file,skip_header=start,max_rows=row,usecols=(1,2,3),unpack=True)
    r = []
    a_count = int(len(rx))
    for c in range(a_count):
        r.append([])
    for i in range(a_count):
        r[i].append(rx[i])
        r[i].append(ry[i])
        r[i].append(rz[i])
    return r

def pbc_dist(r1,r2,L):
    drij = r1 - r2
    drij = drij - L*np.round(drij/L)
    dr = np.sqrt(drij[0]**2. + drij[1]**2. + drij[2]**2.)
    return dr

def lb_mix(e1,e2,s1,s2):
    e = np.sqrt(e1*e2)
    s = (s1+s2)/2
    return e,s

def calc_LJ(dr, epsilon, sigma):
    dsr6 = (sigma/dr)**6.
    dsr12 = dsr6**2
    return 4*epsilon*(dsr12-dsr6)

def calc_COUL(dr, q1, q2):
    C = 332.0739
    return C*q1*q2/dr

def calc_totE(mol,typ,Q,r,L,eps,sig,cutoff):
    e_LJ = 0.0
    e_COUL = 0.0
    for mol1 in range(len(mol)-1):
        for a_m1 in range(len(mol[mol1])):
            n1 = mol[mol1][a_m1]
            q1 = Q[n1]
            for mol2 in range(mol1+1,len(mol)):
                for a_m2 in range(len(mol[mol2])):
                    n2 = mol[mol2][a_m2]
                    q2 = Q[n2]
                    dr = pbc_dist(np.asarray(r[n1]),np.asarray(r[n2]),np.asarray(L))
                    e,s = lb_mix(eps[typ[n1]-1],eps[typ[n2]-1],
                                 sig[typ[n1]-1],sig[typ[n2]-1])
                    if dr < cutoff:
                        LJ = calc_LJ(dr,e,s)
                        e_LJ = e_LJ + LJ
                        COUL = calc_COUL(dr,q1,q2)
                        e_COUL = e_COUL + COUL
    return e_LJ, e_COUL

def write_E():
    return

if __name__ == "__main__":
    import sys
    import numpy as np
    
    # to run type "calc_init_energy.py in.*** data.*** traj.xyz(1 run)" 
    # still working on the cutoff- adjust for now
    inpfile = str(sys.argv[1])
    datafile = str(sys.argv[2])    
    trajfile = str(sys.argv[3])    

    eps, sig = read_paircoeffs(inpfile)
    mol, typ, Q, L = read_data(datafile)
    r = read_frame(trajfile)
    cutoff = 2.5
    e_LJ, e_COUL = calc_totE(mol,typ,Q,r,L,eps,sig,cutoff)
    print(e_LJ, e_COUL)
