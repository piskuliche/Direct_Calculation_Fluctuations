#!/usr/bin/env python
""" 
I realize this isn't a fortran program - but it was much simpler to write.
"""
import numpy as np
import time


# Input Parameters

# SPC/E Water Model

qh = 0.42380 # e
qo = -0.84760 # e
eps = 0.15535 # Kcal/mol
sig = 3.16600 # Angstroms
doh = 1.0 # Angstroms

def LJ(dr,eps,sig):
    """
    This function calculates the lennard jones potential energy for a particular epsilon and sigma
    """
    r6 = (sig/dr)**6.
    Ulj = 4*eps*(r6**2.-r6)
    return Ulj

def Coul(dr,q1,q2):
    """
    This function calculates the real-space electrostatic energy
    """
    C=332.06371
    Ucoul = C*q1*q2/dr
    return Ucoul

def calc_lj(mol1,mol2,eps,sig):
    """
    This calls the calc lj script
    """
    dr = min_dist(rO[mol1],rO[mol2])[0]
    return LJ(dr,eps,sig)

def calc_coul(mol1,mol2,qo,qh):
    """
    This calls the electrostatic calculation script
    """
    terms = []
    # Note - there are more terms because there are many sites which have interactions! 
    # 2 water molecules, 3 atoms each - 9 total terms
    # The donor and acceptor are tagged                 D  A
    terms.append([min_dist(rO[mol1],rO[mol2])[0],qo,qo]) # Ox Ox
    terms.append([min_dist(rO[mol1],r1[mol2])[0],qo,qh]) # Ox H1
    terms.append([min_dist(rO[mol1],r2[mol2])[0],qo,qh]) # Ox H2
    terms.append([min_dist(r1[mol1],rO[mol2])[0],qh,qo]) # H1 Ox
    terms.append([min_dist(r2[mol1],rO[mol2])[0],qh,qo]) # H2 Ox
    terms.append([min_dist(r1[mol1],r1[mol2])[0],qh,qh]) # H1 H1
    terms.append([min_dist(r2[mol1],r1[mol2])[0],qh,qh]) # H2 H1
    terms.append([min_dist(r1[mol1],r2[mol2])[0],qh,qh]) # H1 H2
    terms.append([min_dist(r2[mol1],r2[mol2])[0],qh,qh]) # H2 H2
    Ucoul=0
    # Sums up the terms
    for term in terms:
        Ucoul += Coul(term[0],term[1],term[2])
    return Ucoul


def min_dist(a, b):
    """
    This calculates the minimum image distance between two 3D coordinates, stored 
    in a and b. This outputs the overall final distance as a float.
    """
    dr      = []
    shift   = []
    dsq     = 0.0
    for i in range(3):
        shift.append(L*round((a[i]-b[i])/L))
        dr.append(a[i]-b[i]-shift[i])
        dsq += dr[i]**2.
    drfinal = np.sqrt(dsq)
    return drfinal, dr

def calc_ang(dOO,dHO):
    """
    Calculates the hbond jump angle.
    Note - doh is described aboved
    """
    dHOO = np.arccos(np.around((doh**2. + dOO**2. - dHO**2.)/(2.*doh*dOO),4))
    return dHOO
    
def read_frames(filename):
    """
    Reads the frames within filename.
    """
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

def is_it_hbonded(mol1, mol2, rO, r1, r2, t, orig):
    """
    Outputs whether it is hbonded or not.
    """
    hbonded = 0 # default not hbonded
    ehat = [0.0, 0.0, 0.0] # default not hbonded
    # Calculates OO Distance
    dOO, droo = min_dist(rO[mol1+t],rO[mol2+t])
    if dOO == 0.0: print(dOO,mol1,mol2)
    if dOO < rOO_max:
        # Calculates the distance of each OH on molecule 1 from the Oxygen
        dHO1 = min_dist(r1[mol1+t], rO[mol2+t])[0]
        dHO2 = min_dist(r2[mol1+t], rO[mol2+t])[0]
        if orig == 1 or orig == -1:
            if dHO1 < rHO_max:
                # Calculates the Angle of HOO from the donor and acceptor (first OH)
                dHOO1 = np.degrees(calc_ang(dOO, dHO1))
                if dHOO1 < ang_max:
                    hbonded = 1
        if orig == 2 or orig == -1:
            if dHO2 < rHO_max:
                # Calculates the Angle of HOO from the donor and acceptor (second OH)
                dHOO2 = np.degrees(calc_ang(dOO, dHO2))
                if dHOO2 < ang_max:
                    hbonded = 2
    ehat = np.divide(droo,dOO)
    return hbonded, ehat

def calc_cn(e_o, e_t):
    edote = np.dot(e_o, e_t)
    c1 = edote
    c2 = 0.5*(3*edote**2. - 1)
    c3 = 0.5*(5*edote**3-3*edote)
    return c1, c2, c3

def calc_jumpang(OH,rO,r1,r2,timeindex):
    """
    Calculates the jump angle after a switch
    """
    mol1, molo, mol2 = OH[0],OH[1],OH[6]
    dro_old, oldvec = min_dist(rO[mol1+timeindex],rO[molo+timeindex])
    dro_new, newvec = min_dist(rO[mol1+timeindex],rO[mol2+timeindex])
    e_old = np.divide(oldvec, dro_old)
    e_new = np.divide(newvec, dro_new)
    theta = np.arccos(np.dot(e_old, e_new))
    return theta

tstart =time.time()

# Read Corr Calc Input File
inpfile='corr_calc.in'
inp = open(inpfile, 'r')
lines = []
for line in inp:
    lines.append(line)
nfile  = str(lines[1]).strip()
ntimes = int(lines[3].split()[0])
dt     = float(lines[3].split()[1])
volume = float(lines[5])
mol_name = str(lines[7]).strip()
inp.close()

L = volume**(1/3.)

# Definition of the HBond Criteria
rOO_max = 3.1
rHO_max = 2.0
ang_max = 20.

# Trajectory File Name
filename = 'traj_'+str(nfile)+'_'+str(mol_name)+'.xyz'

# Read the frames
print('Reading frames')
nmols, rO, r1, r2 = read_frames(filename)
# Determine initial HBONDS (time 0)

"""
Form of the OHs vector:
    This has multiple elements that are important
    The first index provides the row of interest, corresponding to a particular OH
    The second index provides the column of interest which is separated into a few different parts
    moloh mola,old on? moldonor1 moldonor2 mola,new
"""
OHs,ehat_init = [],[]

print('Calculating initial hbonds')
for mol1 in range(nmols):
    for mol2 in range(nmols):
        if mol1 != mol2:
            hbnd, ehat = is_it_hbonded(mol1, mol2, rO, r1, r2,0,-1) # hbnd has options 0 (not hbonded) 1 (first oh hbonded) and 2 (second oh hbonded)
            if hbnd != 0:
                OHs.append([mol1,mol2,hbnd])
                ehat_init.append(ehat)

# This section checks for donors (up to a maximum of 3)
for OH1 in range(len(OHs)):
    count = 0
    if count < 2: # Maximum of 3 donors planned
        for OH2 in range(len(OHs)):
            if OHs[OH2][0] == OHs[OH1][1]:
                OHs[OH2].append(OHs[OH1][0])
                count+=1
# This section fills out the donors with -1's where there are no donors
for OH1 in range(len(OHs)):
    if len(OHs[OH1]) == 3:
        OHs[OH1].append(-1)
        OHs[OH1].append(-1)
        OHs[OH1].append(-1)
    elif len(OHs[OH1]) == 4:
        OHs[OH1].append(-1)
        OHs[OH1].append(-1)
    elif len(OHs[OH1]) == 5:
        OHs[OH1].append(-1)
    


print("There are %s OHs" % len(OHs))

# Sets frequency to print progress.
nval = int(ntimes/10.)

crp = np.zeros((len(OHs),ntimes))
c1,c2,c3  = np.zeros((len(OHs), ntimes)),np.zeros((len(OHs), ntimes)),np.zeros((len(OHs), ntimes))
norm_t = np.zeros(ntimes)
theta = []
steps = [0]
m = open('normoldcode.dat','w')
# Loop over time, checks if hbonded still
for n in range(1,ntimes):
    if n%nval == 0:
        print("Step Reached: %d" % n)
    steps.append(n)
    timeindex = n*nmols
    for OH in range(len(OHs)):
        crp[OH][n] = crp[OH][n-1]
        mol1 = OHs[OH][0]
        for mol2 in range(nmols):
            hbnd = 0
            if mol1 != mol2 and OHs[OH][2] != 0:
                hbnd, ehat = is_it_hbonded(mol1, mol2, rO, r1, r2, timeindex,OHs[OH][2])
                if hbnd != 0:
                    if OHs[OH][1] == mol2: # Same Acceptor
                        crp[OH][n]+=0
                        # Calculates the oo cn
                        c1[OH][n],c2[OH][n],c3[OH][n]=calc_cn(ehat_init[OH],ehat)
                        # Time dependent normalization
                        norm_t[n] += 1.0
                    else: # New Acceptor
                        crp[OH][n]+=1
                        c1[mol1][n]=0.0
                        c2[mol1][n]=0.0
                        c3[mol1][n]=0.0
                        OHs[OH][2]=0
                        OHs[OH].append(mol2)
                        theta.append(calc_jumpang(OHs[OH],rO,r1,r2,timeindex))
                if hbnd == 0 and mol2 == OHs[OH][1]:
                    c1[OH][n],c2[OH][n],c3[OH][n]=calc_cn(ehat_init[OH],ehat)
                    norm_t[n] += 1.0 
    m.write("%s %s\n" % (n,norm_t[n]))
m.close()
# This handles the case where there is no new acceptor 
for OH in range(len(OHs)):
    if len(OHs[OH]) != 7:
        OHs[OH].append(-1) 


print("Jiggying up final calculations")

C1 = np.divide(np.sum(c1,axis=0),norm_t,out=np.zeros_like(np.sum(c1,axis=0)),where=norm_t!=0)
C2 = np.divide(np.sum(c2,axis=0),norm_t,out=np.zeros_like(np.sum(c2,axis=0)),where=norm_t!=0)
C3 = np.divide(np.sum(c3,axis=0),norm_t,out=np.zeros_like(np.sum(c3,axis=0)),where=norm_t!=0)
C1[0],C2[0],C3[0]=1.0,1.0,1.0

CRP = np.average(crp,axis=0)

UCRP={"LJAold":np.zeros((len(OHs),ntimes)), "LJAnew":np.zeros((len(OHs),ntimes)), "LJD":np.zeros((len(OHs),ntimes)),"CAold":np.zeros((len(OHs),ntimes)),"CAnew":np.zeros((len(OHs),ntimes)), "CD":np.zeros((len(OHs),ntimes))}
UC2={"LJAold":np.zeros((len(OHs),ntimes)), "LJAnew":np.zeros((len(OHs),ntimes)), "LJD":np.zeros((len(OHs),ntimes)),"CAold":np.zeros((len(OHs),ntimes)),"CAnew":np.zeros((len(OHs),ntimes)), "CD":np.zeros((len(OHs),ntimes))}

Ulj={"Aold":[],"Anew":[], "D":[]}
UCoul={"Aold":[],"Anew":[], "D":[]}


for OH in range(len(OHs)):
    target, oldaccept, donor1, donor2, donor3, newaccept = OHs[OH][0],OHs[OH][1],OHs[OH][3],OHs[OH][4],OHs[OH][5], OHs[OH][6]
    # Old Acceptor
    Ulj["Aold"].append(calc_lj(target,oldaccept,eps,sig))
    UCRP["LJAold"][OH]=np.multiply(Ulj["Aold"][OH],crp[OH])
    UC2["LJAold"][OH]=np.multiply(Ulj["Aold"][OH],c2[OH])
    UCoul["Aold"].append(calc_coul(target,oldaccept,qo,qh))
    UCRP["CAold"][OH]=np.multiply(UCoul["Aold"][OH],crp[OH])
    UC2["CAold"][OH]=np.multiply(UCoul["Aold"][OH],c2[OH])
    # New Acceptor
    if newaccept != -1:
        Ulj["Anew"].append(calc_lj(target,newaccept,eps,sig))
        UCoul["Anew"].append(calc_coul(target,newaccept,qo,qh))
    else:
        Ulj["Anew"].append(0.0)
        UCoul["Anew"].append(0.0)
    UCRP["CAnew"][OH]=np.multiply(UCoul["Anew"][OH],crp[OH])
    UCRP["LJAnew"][OH]=np.multiply(Ulj["Anew"][OH],crp[OH])
    UC2["CAnew"][OH]=np.multiply(UCoul["Anew"][OH],c2[OH])
    UC2["LJAnew"][OH]=np.multiply(Ulj["Anew"][OH],c2[OH])

    tmpljd=0.0
    # Donors
    Ulj["D"].append(0.0)
    UCoul["D"].append(0.0)
    for item in [donor1,donor2,donor3]:
        if item != -1:
            Ulj["D"][OH]   += calc_lj(target,item,eps,sig)
            UCoul["D"][OH] += calc_coul(target,donor1,qo,qh)
    UCRP["LJD"][OH]=np.multiply(Ulj["D"][OH],crp[OH])
    UCRP["CD"][OH]=np.multiply(UCoul["D"][OH],crp[OH])
    UC2["LJD"][OH]=np.multiply(Ulj["D"][OH],c2[OH])
    UC2["CD"][OH]=np.multiply(UCoul["D"][OH],c2[OH])

avcrp,avc2 = {},{}
for key in UCRP:
    avcrp[key]=np.sum(UCRP[key],axis=0)
    np.savetxt(key+'dcrp_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps,avcrp[key]],fmt="%2.5f")

for key in UC2:
    avc2[key]=np.divide(np.sum(UC2[key],axis=0),norm_t,out=np.zeros_like(np.sum(UC2[key],axis=0)),where=norm_t!=0)
    np.savetxt(key+'dframe_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps,avc2[key]],fmt="%2.5f")


tht,bedge = np.histogram(theta, bins=50,range=(0,3))
split=3/50.0
bins=[]
print(len(tht),len(bins))
for b in range(50):
    bins.append(b*split+split/2)

np.savetxt('theta_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[bins,tht], fmt="%2.5f")
np.savetxt('crp_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps, CRP], fmt="%2.5f")
np.savetxt('hbnd_'+str(nfile)+'_'+str(mol_name)+'.dat',np.c_[steps,norm_t], fmt="%2.5f")
np.savetxt('framec1_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps, C1, norm_t], fmt="%2.5f")
np.savetxt('framec2_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps, C2, norm_t], fmt="%2.5f")
np.savetxt('framec3_'+str(nfile)+'_'+str(mol_name)+'.dat', np.c_[steps, C3, norm_t], fmt="%2.5f")
np.savetxt('LJAold_init.out', np.c_[np.sum(Ulj["Aold"])])
np.savetxt('LJAnew_init.out', np.c_[np.sum(Ulj["Anew"])])
np.savetxt('LJD_init.out', np.c_[np.sum(Ulj["D"])])
np.savetxt('CAold_init.out', np.c_[np.sum(UCoul["Aold"])])
np.savetxt('CAnew_init.out', np.c_[np.sum(UCoul["Anew"])])
np.savetxt('CD_init.out', np.c_[np.sum(UCoul["D"])])

print("time", time.time()-tstart)       
