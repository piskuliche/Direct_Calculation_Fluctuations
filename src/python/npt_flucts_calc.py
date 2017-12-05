# This code is used to calculate the activation energies of diffusion and reorientation based on the fluctuations method.
# It should run for either the NVT or NPT ensembles.

from scipy import stats
import numpy as np
import os
import argparse
from argparse import RawTextHelpFormatter,SUPPRESS

parser = argparse.ArgumentParser(description='''Calculates the fluctuations of NPT Ensemble''',formatter_class=RawTextHelpFormatter)
parser.add_argument('-files', help ="Number Files")
parser.add_argument('-blocks', help ="Number of blocks")
args = parser.parse_args()
nfiles=int(args.files)
nblocks=int(args.blocks)

filesperblock=nfiles/float(nblocks)
t_val=stats.t.ppf(0.95,nblocks-1)/np.sqrt(nblocks)
# Vector initialization
tmsd=[400]
msd=[400]
c2=[400]
dmsd=[400]
ea=[400]
fnames=[nfiles]
ke=[nfiles]
dke=[nfiles]
v=[nfiles]
dv=[nfiles]
e=[nfiles]
de=[nfiles]
lj=[nfiles]
dlj=[nfiles]
coul=[nfiles]
dcoul=[nfiles]
ewald=[nfiles]
dewald=[nfiles]
vol=[nfiles]
dvol=[nfiles]
# msd calcs
dkemsd=[400]
demsd=[400]
dvmsd=[400]
dljmsd=[400]
dcoulmsd=[400]
dewaldmsd=[400]
dvolmsd=[400]
# c1 calcs
dkec1=[400]
dec1=[400]
dvc1=[400]
dljc1=[400]
dcoulc1=[400]
dewaldc1=[400]
dvolc1=[400]
# c2 calcs
dkec2=[400]
dec2=[400]
dvc2=[400]
dljc2=[400]
dcoulc2=[400]
dewaldc2=[400]
dvolc2=[400]
# MSD EA calcs
eake_bl=[[] for x in range(0,nblocks)]
eav_bl=[[] for x in range(0,nblocks)]
eae_bl=[[] for x in range(0,nblocks)]
ealj_bl=[[] for x in range(0,nblocks)]
eacoul_bl=[[] for x in range(0,nblocks)]
eaewald_bl=[[] for x in range(0,nblocks)]
eavol_bl=[[] for x in range(0,nblocks)]
# dmsd initialization
dmsd_ke_bl=[[] for x in range(0,nblocks)]
dmsd_v_bl=[[] for x in range(0,nblocks)]
dmsd_e_bl=[[] for x in range(0,nblocks)]
dmsd_lj_bl=[[] for x in range(0,nblocks)]
dmsd_coul_bl=[[] for x in range(0,nblocks)]
dmsd_ewald_bl=[[] for x in range(0,nblocks)]
dmsd_vol_bl=[[] for x in range(0,nblocks)]


# DC2/DbetaP
c2ke_bl=[[] for x in range(0,nblocks)]
c2v_bl=[[] for x in range(0,nblocks)]
c2e_bl=[[] for x in range(0,nblocks)]
c2lj_bl=[[] for x in range(0,nblocks)]
c2coul_bl=[[] for x in range(0,nblocks)]
c2ewald_bl=[[] for x in range(0,nblocks)]
c2vol_bl=[[] for x in range(0,nblocks)]
# C2
c2_bl=[[] for x in range(0,nblocks)]

for i in range(0,nfiles):
    dke.append(0)
    dv.append(0)
    de.append(0)
    dlj.append(0)
    dcoul.append(0)
    dewald.append(0)
    dvol.append(0)
#read in file names, assign lists of these file names
fnames=np.genfromtxt('file_names', dtype='string', unpack=True)
fmsd=[s+"/msd_"+s+".dat" for s in fnames]
fc1=[s+"/c1_"+s+".dat" for s in fnames]
fc2=[s+"/c2_"+s+".dat" for s in fnames]

#read in energy components (and volume)
ke=np.genfromtxt('ke_init.out',unpack=True)
v=np.genfromtxt('v_init.out', unpack=True)
e=np.genfromtxt('e_init.out', unpack=True)
lj=np.genfromtxt('lj_init.out', unpack=True)
coul=np.genfromtxt('coul_init.out', unpack=True)
ewald=np.genfromtxt('ewald_init.out', unpack=True)
vol=np.genfromtxt('vol_init.out', unpack=True)

# Do block calculation
for block in range(0,nblocks):
    start=int(block*filesperblock)
    end=int((block+1)*filesperblock)
    print start, end
    keavg_bl=np.average(ke[start:end])
    vavg_bl=np.average(v[start:end])
    eavg_bl=np.average(e[start:end])
    ljavg_bl=np.average(lj[start:end])
    coulavg_bl=np.average(coul[start:end])
    ewaldavg_bl=np.average(ewald[start:end])
    volavg_bl=np.average(vol[start:end])
    for j in range(start,end):
        dke[j]=ke[j]-keavg_bl
        dv[j]=v[j]-vavg_bl
        de[j]=e[j]-eavg_bl
        dlj[j]=lj[j]-ljavg_bl
        dcoul[j]=coul[j]-coulavg_bl
        dewald[j]=ewald[j]-ewaldavg_bl
        dvol[j]=vol[j]-volavg_bl
    for j in range(0,len(msd)):
        msd[j]=0
        dkemsd[j]=0
        dvmsd[j]=0
        demsd[j]=0
        dljmsd[j]=0
        dcoulmsd[j]=0
        dewaldmsd[j]=0
        dvolmsd[j]=0 
        c2[j]=0
        dkec2[j]=0
        dvc2[j]=0
        dec2[j]=0
        dljc2[j]=0
        dcoulc2[j]=0
        dewaldc2[j]=0
        dvolc2[j]=0
    for j in range(start,end):
        tmsd=np.genfromtxt(fmsd[j], usecols=1, unpack=True)
        tc2=np.genfromtxt(fc2[j], usecols=1, unpack=True)
        # Multiply by d_fluct MSD
        msd=msd+tmsd
        dkemsd=dkemsd+tmsd*dke[int(j)]
        dvmsd=dvmsd+tmsd*dv[int(j)]
        demsd=demsd+tmsd*de[int(j)]
        dljmsd=dljmsd+tmsd*dlj[int(j)]
        dcoulmsd=dcoulmsd+tmsd*dcoul[int(j)]
        dewaldmsd=dewaldmsd+tmsd*dewald[int(j)]
        dvolmsd=dvolmsd+tmsd*dvol[int(j)]
        # Multiply by d_fluct C2
        c2=c2+tc2
        dkec2=dkec2+tc2*dke[int(j)]
        dvc2=dvc2+tc2*dv[int(j)]
        dec2=dec2+tc2*de[int(j)]
        dljc2=dljc2+tc2*dlj[int(j)]
        dcoulc2=dcoulc2+tc2*dcoul[int(j)]
        dewaldc2=dewaldc2+tc2*dewald[int(j)]
        dvolc2=dvolc2+tc2*dvol[int(j)]
    # Normalize the lists
    # MSD
    msd[:] = [x / float(filesperblock) for x in msd]
    dkemsd[:] = [x / float(filesperblock) for x in dkemsd]
    dvmsd[:] = [x / float(filesperblock) for x in dvmsd]
    demsd[:] = [x / float(filesperblock) for x in demsd]
    dljmsd[:] = [x / float(filesperblock) for x in dljmsd]
    dcoulmsd[:] = [x / float(filesperblock) for x in dcoulmsd]
    dewaldmsd[:] = [x / float(filesperblock) for x in dewaldmsd]
    dvolmsd[:] = [x / float(filesperblock) for x in dvolmsd]
    # C2
    c2[:] = [x / float(filesperblock) for x in c2]
    dkec2[:] = [x / float(filesperblock) for x in dkec2]
    dvc2[:] = [x / float(filesperblock) for x in dvc2]
    dec2[:] = [x / float(filesperblock) for x in dec2]
    dljc2[:] = [x / float(filesperblock) for x in dljc2]
    dcoulc2[:] = [x / float(filesperblock) for x in dcoulc2]
    dewaldc2[:] = [x / float(filesperblock) for x in dewaldc2]
    dvolc2[:] = [x / float(filesperblock) for x in dvolc2]
    # Activation Calculations MSD
    eake=dkemsd[1:]/msd[1:]
    eav=dvmsd[1:]/msd[1:]
    eae=demsd[1:]/msd[1:]
    ealj=dljmsd[1:]/msd[1:]
    eacoul=dcoulmsd[1:]/msd[1:]
    eaewald=dewaldmsd[1:]/msd[1:]
    eavol=dvolmsd[1:]/msd[1:]
    # Saves to sep list
    eake_bl[block]=list(eake)
    eav_bl[block]=list(eav)
    eae_bl[block]=list(ealj)
    eacoul_bl[block]=list(eacoul)
    ealj_bl[block]=list(ealj)
    eaewald_bl[block]=list(eaewald)
    eavol_bl[block]=list(eavol)
    # Saves msd_info to sep list
    dmsd_ke_bl[block]=list(dkemsd[1:])
    dmsd_v_bl[block]=list(dvmsd[1:])
    dmsd_e_bl[block]=list(demsd[1:])
    dmsd_lj_bl[block]=list(dljmsd[1:])
    dmsd_coul_bl[block]=list(dcoulmsd[1:])
    dmsd_ewald_bl[block]=list(dewaldmsd[1:])
    dmsd_vol_bl[block]=list(dvolmsd[1:])
    # Save c2 to sep list
    c2_bl[block]=list(c2)
    c2ke_bl[block]=list(dkec2)
    c2v_bl[block]=list(dvc2)
    c2e_bl[block]=list(dec2)
    c2coul_bl[block]=list(dcoulc2)
    c2lj_bl[block]=list(dljc2)
    c2ewald_bl[block]=list(dewaldc2)
    c2vol_bl[block]=list(dvolc2)
    # print the block values
    np.savetxt('bl_'+str(block)+'_msd.dat',np.c_[eake,eav,eae,ealj,eacoul,eaewald,eavol], fmt='%s')
    np.savetxt('bl_'+str(block)+'_c2.dat',np.c_[dkec2,dvc2,dec2,dljc2,dcoulc2,dewaldc2,dvolc2], fmt='%s')
    # Zero the things of range (start,end)
    for i in range(start,end):
        dke[i]=0
        dv[i]=0
        de[i]=0
        dlj[i]=0
        dcoul[i]=0
        dewald[i]=0
        dvol[i]=0
    for j in range(0,len(msd)):
        msd[j]=0
        dkemsd[j]=0
        dvmsd[j]=0
        demsd[j]=0
        dljmsd[j]=0
        dcoulmsd[j]=0
        dewaldmsd[j]=0
        dvolmsd[j]=0
    for k in range(0,len(eake)):
        eake[k]=0
        eav[k]=0
        eae[k]=0
        ealj[k]=0
        eacoul[k]=0
        eaewald[k]=0
        eavol[k]=0
# Do total calculation
keavg=np.average(ke)
vavg=np.average(v)
eavg=np.average(e)
ljavg=np.average(lj)
coulavg=np.average(coul)
ewaldavg=np.average(ewald)
volavg=np.average(vol)
print keavg, vavg, eavg, ljavg, coulavg, ewaldavg, volavg
# Zeroes the msd Vectors and recalculates the d_value lists
for j in range(0,nfiles):
    dke[j]=ke[j]-keavg
    dv[j]=v[j]-vavg
    de[j]=e[j]-eavg
    dlj[j]=lj[j]-ljavg
    dcoul[j]=coul[j]-coulavg
    dewald[j]=ewald[j]-ewaldavg
    dvol[j]=vol[j]-volavg
for j in range(0,len(msd)):    
    msd[j]=0
    dkemsd[j]=0
    dvmsd[j]=0
    demsd[j]=0
    dljmsd[j]=0
    dcoulmsd[j]=0
    dewaldmsd[j]=0
    dvolmsd[j]=0
    c2[j]=0
    dkec2[j]=0
    dvc2[j]=0
    dec2[j]=0
    dljc2[j]=0
    dcoulc2[j]=0
    dewaldc2[j]=0
    dvolc2[j]=0
for j in range(0,nfiles):
    tmsd=np.genfromtxt(fmsd[j], usecols=1, unpack=True)
    tc2=np.genfromtxt(fc2[j], usecols=1, unpack=True)
    # Multiply MSD by d_fluct
    msd=msd+tmsd
    dkemsd=dkemsd+tmsd*dke[int(j)]
    dvmsd=dvmsd+tmsd*dv[int(j)]
    demsd=demsd+tmsd*de[int(j)]
    dljmsd=dljmsd+tmsd*dlj[int(j)]
    dcoulmsd=dcoulmsd+tmsd*dcoul[int(j)]
    dewaldmsd=dewaldmsd+tmsd*dewald[int(j)]
    dvolmsd=dvolmsd+tmsd*dvol[int(j)]
    # Multiply C2 by d_fluct
    c2=c2+tc2
    dkec2=dkec2+tc2*dke[int(j)]
    dvc2=dvc2+tc2*dv[int(j)]
    dec2=dec2+tc2*de[int(j)]
    dljc2=dljc2+tc2*dlj[int(j)]
    dcoulc2=dcoulc2+tc2*dcoul[int(j)]
    dewaldc2=dewaldc2+tc2*dewald[int(j)]
    dvolc2=dvolc2+tc2*dvol[int(j)]  
# Normalize Results
#myList[:] = [x / myInt for x in myList]
msd[:] = [x / float(nfiles) for x in msd]
dkemsd[:] = [x / float(nfiles) for x in dkemsd]
dvmsd[:] = [x / float(nfiles) for x in dvmsd]
demsd[:] = [x / float(nfiles) for x in demsd]
dljmsd[:] = [x / float(nfiles) for x in dljmsd]
dcoulmsd[:] = [x / float(nfiles) for x in dcoulmsd]
dewaldmsd[:] = [x / float(nfiles) for x in dewaldmsd]
dvolmsd[:] = [x / float(nfiles) for x in dvolmsd]
# C2
c2[:] = [x / float(nfiles) for x in c2]
dkec2[:] = [x / float(nfiles) for x in dkec2]
dvc2[:] = [x / float(nfiles) for x in dvc2]
dec2[:] = [x / float(nfiles) for x in dec2]
dljc2[:] = [x / float(nfiles) for x in dljc2]
dcoulc2[:] = [x / float(nfiles) for x in dcoulc2]
dewaldc2[:] = [x / float(nfiles) for x in dewaldc2]
dvolc2[:] = [x / float(nfiles) for x in dvolc2]
# Activation Calculations
eake=dkemsd[1:]/msd[1:]
eav=dvmsd[1:]/msd[1:]
eae=demsd[1:]/msd[1:]
ealj=dljmsd[1:]/msd[1:]
eacoul=dcoulmsd[1:]/msd[1:]
eaewald=dewaldmsd[1:]/msd[1:]
eavol=dvolmsd[1:]/msd[1:]
# Save Data
np.savetxt('tot_msd_ea.dat', np.c_[eake,eav,eae,ealj,eacoul,eaewald,eavol], fmt='%s')
np.savetxt('tot_dc2dbp.dat', np.c_[dkec2,dvc2,dec2,dljc2,dcoulc2,dewaldc2,dvolc2], fmt='%s')
np.savetxt('tot_msd.dat', msd, fmt=['%.4f'])
np.savetxt('tot_dmsd.dat', np.c_[dkemsd,dvmsd,demsd,dljmsd,dcoulmsd,dewaldmsd,dvolmsd], fmt='%s')

# Save Component Distributions
np.savetxt('dke.dat', dke, fmt=['%.4f'])
np.savetxt('dv.dat', dv, fmt=['%.4f'])
np.savetxt('de.dat', de, fmt=['%.4f'])
np.savetxt('dlj.dat', dlj, fmt=['%.4f'])
np.savetxt('dcoul.dat', dcoul, fmt=['%.4f'])
np.savetxt('dewald.dat', dewald, fmt=['%.4f'])
np.savetxt('dvol.dat', dvol, fmt=['%.4f'])

# calculate uncertainty
    # MSD
err_eake=np.array(eake_bl).std(0)
err_eav=np.array(eav_bl).std(0)
err_eae=np.array(eae_bl).std(0)
err_ealj=np.array(ealj_bl).std(0)
err_eacoul=np.array(eacoul_bl).std(0)
err_eaewald=np.array(eaewald_bl).std(0)
err_eavol=np.array(eavol_bl).std(0)
    #DMSD
err_dkemsd=np.array(dmsd_ke_bl).std(0)
err_dvmsd=np.array(dmsd_v_bl).std(0)
err_demsd=np.array(dmsd_e_bl).std(0)
err_dljmsd=np.array(dmsd_lj_bl).std(0)
err_dcoulmsd=np.array(dmsd_coul_bl).std(0)
err_dewaldmsd=np.array(dmsd_ewald_bl).std(0)
err_dvolmsd=np.array(dmsd_vol_bl).std(0)

    # C2
err_c2ke=np.array(c2ke_bl).std(0)
err_c2v=np.array(c2v_bl).std(0)
err_c2e=np.array(c2e_bl).std(0)
err_c2lj=np.array(c2lj_bl).std(0)
err_c2coul=np.array(c2coul_bl).std(0)
err_c2ewald=np.array(c2ewald_bl).std(0)
err_c2vol=np.array(c2vol_bl).std(0)

# Do confidence interval
    # MSD
err_eake = [ x * t_val for x in err_eake ]
err_eav = [ x * t_val for x in err_eav ]
err_eae = [ x * t_val for x in err_eae ]
err_ealj = [ x * t_val for x in err_ealj ]
err_eacoul = [ x * t_val for x in err_eacoul ]
err_eaewald = [ x * t_val for x in err_eaewald ]
err_eavol = [ x * t_val for x in err_eavol ]
    #DMSD
err_dkemsd = [ x * t_val for x in err_dkemsd ]
err_dvmsd = [ x * t_val for x in err_dvmsd ]
err_demsd = [ x * t_val for x in err_demsd ]
err_dljmsd = [ x * t_val for x in err_dljmsd ]
err_dcoulmsd = [ x * t_val for x in err_dcoulmsd ]
err_dewaldmsd = [ x * t_val for x in err_dewaldmsd ]
err_dvolmsd = [ x * t_val for x in err_dvolmsd ]
    # C2
err_c2ke = [ x * t_val for x in err_c2ke ]
err_c2v = [ x * t_val for x in err_c2v ]
err_c2e = [ x * t_val for x in err_c2e ]
err_c2lj = [ x * t_val for x in err_c2lj ]
err_c2coul = [ x * t_val for x in err_c2coul ]
err_c2ewald = [ x * t_val for x in err_c2ewald ]
err_c2vol = [ x * t_val for x in err_c2vol ] 

# Calculate Time
time=[]
for i in range(0,len(err_eake)+1):
    time.append(i*.05)

np.savetxt('ea_msd_ke.dat', np.c_[time[1:],eake,err_eake], fmt='%s')
np.savetxt('ea_msd_v.dat', np.c_[time[1:],eav,err_eav], fmt='%s')
np.savetxt('ea_msd_e.dat', np.c_[time[1:],eae,err_eae], fmt='%s')
np.savetxt('ea_msd_coul.dat', np.c_[time[1:],eacoul,err_eacoul], fmt='%s')
np.savetxt('ea_msd_lj.dat', np.c_[time[1:],ealj,err_ealj], fmt='%s')
np.savetxt('ea_msd_vol.dat', np.c_[time[1:],eavol,err_eavol], fmt='%s')

np.savetxt('dke_msd_tot.dat', np.c_[time[1:], dkemsd[1:], err_dkemsd], fmt='%s')
np.savetxt('dv_msd_tot.dat', np.c_[time[1:], dvmsd[1:], err_dvmsd], fmt='%s')
np.savetxt('de_msd_tot.dat', np.c_[time[1:], demsd[1:], err_demsd], fmt='%s')
np.savetxt('dlj_msd_tot.dat', np.c_[time[1:], dljmsd[1:], err_dljmsd], fmt='%s')
np.savetxt('dcoul_msd_tot.dat', np.c_[time[1:], dcoulmsd[1:], err_dcoulmsd], fmt='%s')
np.savetxt('dewald_msd_tot.dat', np.c_[time[1:], dewaldmsd[1:], err_dewaldmsd], fmt='%s')
np.savetxt('dvol_msd_tot.dat', np.c_[time[1:], dvolmsd[1:], err_dvolmsd], fmt='%s')





np.savetxt('dc2_dbp_ke.dat', np.c_[time,dkec2, err_c2ke], fmt='%s')
np.savetxt('dc2_dbp_v.dat', np.c_[time,dvc2, err_c2v], fmt='%s')
np.savetxt('dc2_dbp_e.dat', np.c_[time,dec2, err_c2e], fmt='%s')
np.savetxt('dc2_dbp_lj.dat', np.c_[time,dljc2, err_c2lj], fmt='%s')
np.savetxt('dc2_dbp_coul.dat', np.c_[time,dcoulc2, err_c2coul], fmt='%s')
np.savetxt('dc2_dbp_ewald.dat', np.c_[time,dewaldc2, err_c2ewald], fmt='%s')
np.savetxt('dc2_dbp_vol.dat', np.c_[time,dvolc2, err_c2vol], fmt='%s')
# Save stuff for fitting c2
np.savetxt('c2_total_result.dat',np.c_[time,c2],fmt='%s')
np.savetxt('dc2_dbp_e_total_result.dat', np.c_[time,dec2],fmt='%s')
np.savetxt('dc2_dbp_vol_total_result.dat', np.c_[time, dvolc2], fmt='%s')
np.savetxt('time.dat',time,fmt='%s')
for i in range(0,nblocks):
    np.savetxt('bl_'+str(int(i))+'_c2_val.dat', c2_bl[i], fmt='%s')
    np.savetxt('bl_'+str(int(i))+'_dc2_dbp_val_e.dat', c2e_bl[i], fmt='%s')
    np.savetxt('bl_'+str(int(i))+'_dc2_dbp_val_vol.dat', c2vol_bl[i], fmt='%s')

