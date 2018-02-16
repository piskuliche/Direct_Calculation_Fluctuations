# This code is used to calculate the activation energies of diffusion and reorientation based on the fluctuations method.
# It should run for either the NVT or NPT ensembles.

from scipy import stats
from scipy.optimize import curve_fit
from scipy.integrate import quad
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from argparse import RawTextHelpFormatter

def Corr_Fit(t,A1,A2,A3,k1,k2,k3):
    return A1*np.exp(-k1*t)+A2*np.exp(-k2*t)+A3*np.exp(-k3*t)

def dCorr_Fit(t,dA1,dA2,dA3,dk1,dk2,dk3):
    return dA1*np.exp(-k1*t)+dA2*np.exp(-k2*t)+dA3*np.exp(-k3*t)-A1*dk1*t*np.exp(-k1*t)-A2*dk2*t*np.exp(-k2*t)-A3*dk3*t*np.exp(-k3*t)

def d2Corr_Fit(t,d2A1, d2A2, d2A3, d2k1, d2k2, d2k3):
    return (d2A1-t*(2*dA1*dk1+A1*d2k2)+t**2*(A1*dk1*dk1))*np.exp(-k1*t)+(d2A2-t*(2*dA2*dk2+A2*d2k2)+t**2*(A2*dk2*dk2))*np.exp(-k2*t)+(d2A3-t*(2*dA3*dk3+A3*d2k3)+t**2*(A3*dk3*dk3))*np.exp(-k3*t)

def Integrate(x,y):
    return np.trapz(y,x)

def Int2_Prediction(x,T,itau,idc2,id2c2):
    return (1.0/itau)*(1.0+(idc2/itau)*(x-(1.0/(kb*T)))+0.5*(2.0*(idc2/itau)**2.0-id2c2/itau)*(x-(1.0/(kb*T)))**2)

def Error(value):
    err_value = np.array(value).std(0) * t_val
    return err_value

def E_Mult_Prop(a,ea,b,eb):
    return a*b*np.sqrt((ea/a)**2+(eb/b)**2)

# Read in Arguments
parser = argparse.ArgumentParser(description='''Calculates the c2 fit''', formatter_class=RawTextHelpFormatter)
parser.add_argument('-inp', help="Input File")
parser.add_argument('-blocks', help="Number of blocks")
parser.add_argument('-mol', help ="Molecule Name")
parser.add_argument('-temp', help = "Temperature (in Kelvin)")
args = parser.parse_args()

inputfile = str(args.inp)
nblocks = int(args.blocks)
mol_name = str(args.mol)
T=float(args.temp)

#Physical Constants
kb=0.0019872041

# Read in input file
inp_names, inp_cols=np.genfromtxt(inputfile, usecols=(0,1), dtype=(str,int),unpack=True)
print inp_names
print inp_cols

# Calculate important quantities
t_val=stats.t.ppf(0.95,nblocks-1)/np.sqrt(nblocks)

# Initialize Block Arrays
A1_bl=[[] for x in range(0,nblocks)]
A2_bl=[[] for x in range(0,nblocks)]
A3_bl=[[] for x in range(0,nblocks)]
dA1_bl=[[] for x in range(0,nblocks)]
dA2_bl=[[] for x in range(0,nblocks)]
dA3_bl=[[] for x in range(0,nblocks)]
d2A1_bl=[[] for x in range(0,nblocks)]
d2A2_bl=[[] for x in range(0,nblocks)]
d2A3_bl=[[] for x in range(0,nblocks)]
k1_bl=[[] for x in range(0,nblocks)]
k2_bl=[[] for x in range(0,nblocks)]
k3_bl=[[] for x in range(0,nblocks)]
dk1_bl=[[] for x in range(0,nblocks)]
dk2_bl=[[] for x in range(0,nblocks)]
dk3_bl=[[] for x in range(0,nblocks)]
d2k1_bl=[[] for x in range(0,nblocks)]
d2k2_bl=[[] for x in range(0,nblocks)]
d2k3_bl=[[] for x in range(0,nblocks)]
int_tau_bl=[[] for x in range(0,nblocks)]
int_dc2_bl=[[] for x in range(0,nblocks)]
int_d2c2_bl=[[] for x in range(0,nblocks)]


beta=[[] for x in range(0,200)]
order2pred=[[] for x in range(0,200)]
order2prederr=[[] for x in range(0,200)]
order2pred_bl=[[] for x in range(0,nblocks)]

iindex=0
jindex=0

for item1 in inp_names:
    jindex=0
    for item2 in inp_names:
        if jindex >= iindex:
            for block in range(0,nblocks):              
                time, cab, dcab = np.genfromtxt('bl_' + str(block) + '_' + item1 + '_' + str(mol_name) + '_c2.dat', usecols=(0,1,2), unpack=True)
                popt_cab, pcov_cab = curve_fit(Corr_Fit, time, cab, p0=(.2,.3,.5,.05,.2,1))
                A=popt_cab[:3]
                k=popt_cab[3:]
                ksrt=sorted(k)
                A1=A[np.where(k==ksrt[0])]
                A2=A[np.where(k==ksrt[1])]
                A3=A[np.where(k==ksrt[2])]
                k1=ksrt[0]
                k2=ksrt[1]
                k3=ksrt[2]
                popt_dcab, pcov_dcab = curve_fit(dCorr_Fit,time, dcab, p0=(0.2,0.3,0.4,0.2,0.3,0.4))
                dA=popt_dcab[:3]
                dk=popt_dcab[3:]
                dA1=dA[0]
                dA2=dA[1]
                dA3=dA[2]
                dk1=dk[0]
                dk2=dk[1]
                dk3=dk[2]
                # Write single item info to file
                A1_bl[block] = A1
                A2_bl[block] = A2
                A3_bl[block] = A3
                k1_bl[block] = k1
                k2_bl[block] = k2
                k3_bl[block] = k3
                dA1_bl[block] = dA1
                dA2_bl[block] = dA2
                dA3_bl[block] = dA3
                dk1_bl[block] = dk1
                dk2_bl[block] = dk2
                dk3_bl[block] = dk3                               
                d2cab =np.genfromtxt('bl_' + str(block) + '_' + item1 + '_' + item2 + '_' + str(mol_name) + '_c2.dat', usecols=(3), unpack=True)
                popt_d2cab, pcov_d2cab = curve_fit(d2Corr_Fit, time, d2cab, p0=(0.2,0.3,0.4,0.2,0.3,0.4))
                d2A=popt_d2cab[:3]
                d2k=popt_d2cab[3:]
                d2A1=d2A[0]
                d2A2=d2A[1]
                d2A3=d2A[2]
                d2k1=d2k[0]
                d2k2=d2k[1]
                d2k3=d2k[2]
                d2A1_bl[block] = d2A1
                d2A2_bl[block] = d2A2
                d2A3_bl[block] = d2A3
                d2k1_bl[block] = d2k1
                d2k2_bl[block] = d2k2
                d2k3_bl[block] = d2k3
                int_tau_bl[block] = Integrate(time,cab)
                int_dc2_bl[block] = Integrate(time,dcab)
                int_d2c2_bl[block] = Integrate(time,d2cab)
                for b in range(0,200):
                    beta[b]=0.5+float(b)*0.01
                    order2pred[b]=Int2_Prediction(beta[b],T,int_tau_bl[block], int_dc2_bl[block], int_d2c2_bl[block])
                order2pred_bl[block]=order2pred
                
            # I apologize for the lazy coding in the next section - for clarity:
            # Yes I am intentionally writing all the errors out to files just so they
            # can be read in later in the same program during the total calculation 
            err_A1 = Error(A1_bl)
            err_A2 = Error(A2_bl)
            err_A3 = Error(A3_bl)
            err_dA1 = Error(dA1_bl)
            err_dA2 = Error(dA2_bl)
            err_dA3 = Error(dA3_bl)
            err_d2A1 = Error(d2A1_bl)
            err_d2A2 = Error(d2A2_bl)
            err_d2A3 = Error(d2A3_bl)
            err_k1 = Error(k1_bl)
            err_k2 = Error(k2_bl)
            err_k3 = Error(k3_bl)
            err_dk1 = Error(dk1_bl)
            err_dk2 = Error(dk2_bl)
            err_dk3 = Error(dk3_bl)
            err_d2k1 = Error(d2k1_bl)
            err_d2k2 = Error(d2k2_bl)
            err_d2k3 = Error(d2k3_bl)
            err_int_tau = Error(int_tau_bl)
            err_int_dc2 = Error(int_dc2_bl)
            err_int_d2c2 = Error(int_d2c2_bl)
            order2prederr=np.array(order2pred_bl).std(0)
            order2prederr=[x * t_val for x in order2prederr]

            np.savetxt('fiterr_'+item1+'_'+str(mol_name)+'_c2.dat', np.c_[err_A1, err_A2, err_A3, err_k1, err_k2, err_k3, err_dA1, err_dA2, err_dA3, err_dk1, err_dk2, err_dk3], fmt='%s')
            np.savetxt('fiterr_'+item1+'_'+item2+'_'+str(mol_name)+'_c2.dat', np.c_[err_d2A1, err_d2A2, err_d2A3, err_d2k1, err_d2k2, err_d2k3], fmt='%s')
            np.savetxt('fiterr_'+item1+'_'+item2+'_'+str(mol_name)+'_intc2.dat', np.c_[err_int_tau, err_int_dc2, err_int_d2c2], fmt='%s')
            np.savetxt('fiterr_'+item1+'_'+item2+'_'+str(mol_name)+'_errc2.dat', np.c_[order2prederr], fmt='%s')
        jindex+=1
    iindex+=1

print "Beginning total calculation"
# Begin total Calculation
iindex=0
jindex=0
for item1 in inp_names:
    jindex=0
    for item2 in inp_names:
        if jindex >= iindex:
            time, cab, dcab = np.genfromtxt(item1 + '_' + str(mol_name) + '_c2.dat', usecols=(0,1,3), unpack=True)
            popt_cab, pcov_cab = curve_fit(Corr_Fit, time, cab, p0=(.2,.3,.5,.05,.2,1))
            A=popt_cab[:3]
            k=popt_cab[3:]
            ksrt=sorted(k)
            A1=A[np.where(k==ksrt[0])]
            A2=A[np.where(k==ksrt[1])]
            A3=A[np.where(k==ksrt[2])]
            k1=ksrt[0]
            k2=ksrt[1]
            k3=ksrt[2]
            popt_dcab, pcov_dcab = curve_fit(dCorr_Fit,time, dcab, p0=(0.2,0.3,0.4,0.2,0.3,0.4))
            dA=popt_dcab[:3]
            dk=popt_dcab[3:]
            dA1=dA[0]
            dA2=dA[1]
            dA3=dA[2]
            dk1=dk[0]
            dk2=dk[1]
            dk3=dk[2]
            d2cab =np.genfromtxt(item1 + '_' + item2 + '_' + str(mol_name) + '_c2.dat', usecols=(5), unpack=True)
            popt_d2cab, pcov_d2cab = curve_fit(d2Corr_Fit, time, d2cab, p0=(0.2,0.3,0.4,0.2,0.3,0.4))
            d2A=popt_d2cab[:3]
            d2k=popt_d2cab[3:]
            d2A1=d2A[0]
            d2A2=d2A[1]
            d2A3=d2A[2]
            d2k1=d2k[0]
            d2k2=d2k[1]
            d2k3=d2k[2]

            err_A1, err_A2, err_A3, err_k1, err_k2, err_k3, err_dA1, err_dA2, err_dA3, err_dk1, err_dk2, err_dk3 = np.genfromtxt('fiterr_'+item1+'_'+str(mol_name)+'_c2.dat', usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True)
            err_d2A1, err_d2A2, err_d2A3, err_d2k1, err_d2k2, err_d2k3 = np.genfromtxt('fiterr_'+item1+'_'+item2+'_'+str(mol_name)+'_c2.dat', usecols=(0,1,2,3,4,5), unpack=True)
            err_int_tau, err_int_dc2, err_int_d2c2 = np.genfromtxt('fiterr_'+item1+'_'+item2+'_'+str(mol_name)+'_intc2.dat', usecols=(0,1,2), unpack=True)
            order2prederr = np.genfromtxt('fiterr_'+item1+'_'+item2+'_'+str(mol_name)+'_errc2.dat', usecols=(0), unpack=True)

            int_tau=Integrate(time,cab)
            int_dc2=Integrate(time,dcab)
            int_d2c2=Integrate(time,d2cab)
            for b in range(0,200):
                beta[b]=0.5+float(b)*0.01
                order2pred[b]=Int2_Prediction(beta[b],T,int_tau, int_dc2, int_d2c2)

            filepath=str('fit_results_'+item1+'_'+item2+'_'+str(mol_name)+'_c2.dat')
            fout=open(filepath,'w')
            fout.write("Correlation Function: C2(t)\n")
            fout.write("Item 1: %s Item 2: %s\n" % (item1, item2))
            fout.write("  A1: %s   err: %s\n" % (A1[0],err_A1))
            fout.write("  A2: %s   err: %s\n" % (A2[0], err_A2))
            fout.write("  A3: %s   err: %s\n" % (A3[0], err_A3))
            fout.write("  k1: %s   err: %s\n" % (k1, err_k1))
            fout.write("  k2: %s   err: %s\n" % (k2, err_k2))
            fout.write("  k3: %s   err: %s\n" % (k3, err_k3))
            fout.write(" dA1: %s   err: %s\n" % (dA1, err_dA1))
            fout.write(" dA2: %s   err: %s\n" % (dA2, err_dA2))
            fout.write(" dA3: %s   err: %s\n" % (dA3, err_dA3))
            fout.write(" dk1: %s   err: %s\n" % (dk1, err_dk1))
            fout.write(" dk2: %s   err: %s\n" % (dk2, err_dk2))
            fout.write(" dk3: %s   err: %s\n" % (dk3, err_dk3))
            fout.write("d2A1: %s   err: %s\n" % (d2A1, err_d2A1))
            fout.write("d2A2: %s   err: %s\n" % (d2A2, err_d2A2))
            fout.write("d2A3: %s   err: %s\n" % (d2A3, err_d2A3))
            fout.write("d2k1: %s   err: %s\n" % (d2k1, err_d2k1))
            fout.write("d2k2: %s   err: %s\n" % (d2k2, err_d2k2))
            fout.write("d2k3: %s   err: %s\n" % (d2k3, err_d2k3))
            fout.write("\n")
            fout.write(" Ea Tau 2: %s   err: %s\n" % (dk1*(1.0/k1), E_Mult_Prop(dk1,err_dk1, k1, err_k1)))
            fout.write(" Tau 2 Time: %s +/- %s ps\n" % ((1.0/k1),(err_k1/k1*(1.0/k1))))
            fout.write(" Ea Librational: %s   err: %s\n" % (dk2*(1.0/k2), E_Mult_Prop(dk2,err_dk2, k2, err_k2)))
            fout.write(" Librational Time: %s +/- %s ps\n" % ((1.0/k2),(err_k2/k2*(1.0/k2)))) 
            fout.write(" Ea Inertial: %s   err: %s\n" % (dk3*(1.0/k3), E_Mult_Prop(dk3,err_dk3, k3, err_k3)))
            fout.write(" Inertial Time: %s +/- %s ps\n" % ((1.0/k3), (err_k3/k3*(1.0/err_k3))))
            fout.write("Integrated Times\n")
            fout.write("<tau2> = %s +/- %s\n" % (int_tau,err_int_tau))
            fout.write("intdc2 = %s +/- %s\n" % (int_dc2,err_int_dc2))
            fout.write("intd2c2 = %s +/- %s\n" % (int_d2c2,err_int_d2c2))
            fout.write("End Correlation function\n")
            fout.close()
            cabfit=Corr_Fit(time,*popt_cab)
            dcabfit=dCorr_Fit(time, *popt_dcab)
            d2cabfit=d2Corr_Fit(time, *popt_d2cab)
            print len(order2prederr)
            np.savetxt('total_fit_'+str(item1)+'_'+str(item2)+'_'+str(mol_name)+'_c2.dat', np.c_[time,cabfit, dcabfit,d2cabfit], fmt='%s')
            np.savetxt('int_fit_'+str(item1)+'_'+str(item2)+'_'+str(mol_name)+'_c2.dat', np.c_[beta, order2pred, order2prederr], fmt='%s')
        jindex+=1
    iindex+=1
