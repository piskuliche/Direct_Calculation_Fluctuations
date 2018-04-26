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
    return A1*np.exp(-k1*t)+A2*np.exp(-k2*t)+(1-A1-A2)*np.exp(-k3*t)

def Alt_CorrFit(t,aA1,aA2,aA3,ak1,ak2,ak3, Y):
    return aA1*np.exp(-(t*ak1)**Y)+aA2*np.exp(-t*ak2)+(1-aA1-aA2)*np.exp(-(t*ak3)**2)

def dCorr_Fit(t,dA1,dA2,dA3,dk1,dk2,dk3):
    return dA1*np.exp(-k1*t)+dA2*np.exp(-k2*t)+(0.0-dA1-dA2)*np.exp(-k3*t)-A1*dk1*t*np.exp(-k1*t)-A2*dk2*t*np.exp(-k2*t)-A3*dk3*t*np.exp(-k3*t)

def Alt_dCorrFit(t,adA1,adA2,adA3,adk1,adk2,adk3):
    return (adA1-Y*aA1*t**Y*ak1**(Y-1)*adk1)*np.exp(-(t*ak1)**Y)+(adA2-aA2*adk2*t)*np.exp(-t*ak2)+(adA3-2*t**2*aA3*ak3*adk3)*np.exp(-(t*ak3)**2)

def d2Corr_Fit(t,d2A1, d2A2, d2A3, d2k1, d2k2, d2k3):
    return (d2A1-t*(2*dA1*dk1+A1*d2k1)+t**2*(A1*dk1*dk1))*np.exp(-k1*t)+(d2A2-t*(2*dA2*dk2+A2*d2k2)+t**2*(A2*dk2*dk2))*np.exp(-k2*t)+(d2A3-t*(2*dA3*dk3+A3*d2k3)+t**2*(A3*dk3*dk3))*np.exp(-k3*t)

def d3Corr_Fit(t, d3A1, d3A2, d3A3, d3k1, d3k2, d3k3):

    third_fit_val = 0
    A=[A1,A2,A3]
    k=[k1, k2, k3]
    dA=[dA1, dA2, dA3]
    dk=[dk1, dk2, dk3]
    d2A=[d2A1, d2A2, d2A3]
    d2k=[d2k1, d2k2, d2k3]
    d3A=[d3A1, d3A2, d3A3]
    d3k=[d3k1, d3k2, d3k3]
    for i in range(3):
        third_fit_val += (d3A[i] - t * (3*d2A[i]*dk[i] +3*dA[i]*d2k[i]+A[i]*d3k[i]) + t**2 * (3*A[i]*dk[i]*d2k[i]+3*dA[i]*(dk[i])**2) - t**3 * (A[i] * dk[i]**3))*np.exp(-k[i]*t)
    return third_fit_val

def Integrate(x,y):
    return np.trapz(y,x)

def Order2_Prediction(x,T,itau,idc2,id2c2):
    bzero=(1.0/(kb*T))
    invtau=(1.0/itau)
    return invtau*(1.0+invtau*(idc2)*(x-bzero)+0.5*invtau*(2.0*invtau*(idc2)**2.0-id2c2)*(x-bzero)**2)

def Order3_Prediction(x,T,itau,idc2,id2c2,id3c2):
    bzero=(1.0/(kb*T))
    invtau=(1.0/itau)
    return invtau*(1.0 + invtau*(idc2)*(x-bzero) + 0.5*invtau*(2.0*invtau*(idc2)**2.0 - id2c2)*(x-bzero)**2.0 - (1.0/6.0)*invtau*(6.0*invtau**2*(idc2)**3.0 - 6.0*invtau*idc2*id2c2 + id3c2)*(x-bzero)**3)

def Error(value):
    err_value = np.array(value).std(0) * t_val
    return err_value

def E_Mult_Prop(a,ea,b,eb):
    return a*b*np.sqrt((ea/a)**2+(eb/b)**2)

def E_Divi_Prop(a,ea,b,eb):
    return a/b*np.sqrt((ea/a)**2+(eb/b)**2)

def Ea_Breakdown(A, dA, tau, inttau, Ea):
    comp1=dA*tau/inttau
    comp2=A*tau/inttau*Ea
    return comp1, comp2

# Read in Arguments
parser = argparse.ArgumentParser(description='''Calculates the c2 fit''', formatter_class=RawTextHelpFormatter)
parser.add_argument('-inp', help="Input File")
parser.add_argument('-blocks', help="Number of blocks")
parser.add_argument('-mol', help ="Molecule Name")
parser.add_argument('-temp', help = "Temperature (in Kelvin)")
parser.add_argument('-c2cut', help = "how many points to include in initial cut")
args = parser.parse_args()

inputfile = str(args.inp)
nblocks = int(args.blocks)
mol_name = str(args.mol)
T=float(args.temp)
cut=int(args.c2cut)

#Physical Constants
kb=0.0019872041

# Read in input file
inp_names, inp_cols=np.genfromtxt(inputfile, usecols=(0,1), dtype=(str,int),unpack=True)
print inp_names
print inp_cols

# Calculate important quantities
t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)

sig_arr=np.genfromtxt('real_time.dat', usecols=1, unpack=True)
print len(sig_arr)
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
d3A1_bl=[[] for x in range(0,nblocks)]
d3A2_bl=[[] for x in range(0,nblocks)]
d3A3_bl=[[] for x in range(0,nblocks)]
k1_bl=[[] for x in range(0,nblocks)]
k2_bl=[[] for x in range(0,nblocks)]
k3_bl=[[] for x in range(0,nblocks)]
dk1_bl=[[] for x in range(0,nblocks)]
dk2_bl=[[] for x in range(0,nblocks)]
dk3_bl=[[] for x in range(0,nblocks)]
d2k1_bl=[[] for x in range(0,nblocks)]
d2k2_bl=[[] for x in range(0,nblocks)]
d2k3_bl=[[] for x in range(0,nblocks)]
d3k1_bl=[[] for x in range(0,nblocks)]
d3k2_bl=[[] for x in range(0,nblocks)]
d3k3_bl=[[] for x in range(0,nblocks)]
int_tau_bl=[[] for x in range(0,nblocks)]
int_dc2_bl=[[] for x in range(0,nblocks)]
int_d2c2_bl=[[] for x in range(0,nblocks)]
int_d3c2_bl=[[] for x in range(0,nblocks)]
int_ea_bl=[[] for x in range(0,nblocks)]
tau2_comp1_bl =[[] for x in range(0,nblocks)]
tau2_comp2_bl =[[] for x in range(0,nblocks)]
taul_comp1_bl =[[] for x in range(0,nblocks)]
taul_comp2_bl =[[] for x in range(0,nblocks)]
taui_comp1_bl =[[] for x in range(0,nblocks)]
taui_comp2_bl =[[] for x in range(0,nblocks)]
comp_sum_bl = [[] for x in range(0,nblocks)]

beta=[[] for x in range(0,200)]
order2pred=[[] for x in range(0,200)]
order2prederr=[[] for x in range(0,200)]
order2pred_bl=[[] for x in range(0,nblocks)]
order3pred=[[] for x in range(0,200)]
order3prederr=[[] for x in range(0,200)]
order3pred_bl=[[] for x in range(0,nblocks)]

iindex=0
jindex=0

for item1 in inp_names:
    jindex=0
    for item2 in inp_names:
        if jindex >= iindex:
            item3=item1
            for block in range(0,nblocks):
                # Read in correlation function, 1st derivative
                time, cab, dcab = np.genfromtxt('bl_' + str(block) + '_' + item1 + '_' + str(mol_name) + '_c2.dat', usecols=(0,1,2), unpack=True)
                # Fit correlation function
                popt_cab, pcov_cab = curve_fit(Corr_Fit, time[:cut], cab[:cut], p0=(.2,.5,.5,.002,.3,3),bounds=((0, 0, 0, 0, 0, 0), (1, 1, 1, np.inf, np.inf, np.inf)))
                # Set parameters to variables
                A=popt_cab[:3]
                k=popt_cab[3:]
                ksrt=sorted(k)
                A1=A[np.where(k==ksrt[0])]
                A2=A[np.where(k==ksrt[1])]
                A3=1-A1-A2
                k1=ksrt[0]
                k2=ksrt[1]
                k3=ksrt[2]
                # Fit second correlation function
                print("Block %s First Derivative:" % block)
                popt_dcab, pcov_dcab = curve_fit(dCorr_Fit,time, dcab, p0=(0.2,0.3,0.4,1.0,2.0,3.0))
                # Set parameters to variables
                dA=popt_dcab[:3]
                dk=popt_dcab[3:]
                dA1=dA[0]
                dA2=dA[1]
                dA3=0.0-dA1-dA2
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
                # Read and Fit Second Derivative 
                print("Block %s Second Derivative:" % block)
                d2cab =np.genfromtxt('bl_' + str(block) + '_' + item1 + '_' + item2 + '_' + str(mol_name) + '_c2.dat', usecols=(3), unpack=True)
                popt_d2cab, pcov_d2cab = curve_fit(d2Corr_Fit, time, d2cab, p0=(-1,0.9,3.0,3.3,-20.0,-100.))
                d2A=popt_d2cab[:3]
                d2k=popt_d2cab[3:]
                print d2A
                print d2k
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
                # Read and Fit Third Derivative
                print("Block %s Third Derivative:" % block)
                d3cab = np.genfromtxt('bl_'+str(block)+'_'+item1+'_'+item2+'_'+item3+'_'+mol_name+'_c2.dat', usecols=(1), unpack=True)
                popt_d3cab, pcov_d3cab = curve_fit(d3Corr_Fit, time, d3cab, p0=(0.1,.2,0.3,-0.2,0.3,0.4))
                d3A=popt_d3cab[:3]
                d3k=popt_d3cab[3:]
                d3A1=d3A[0]
                d3A2=d3A[1]
                d3A3=d3A[2]
                d3k1=d3k[0]
                d3k2=d3k[1]
                d3k3=d3k[2]
                d3A1_bl[block] = d3A1
                d3A2_bl[block] = d3A2
                d3A3_bl[block] = d3A3
                d3k1_bl[block] = d3k1
                d3k2_bl[block] = d3k2
                d3k3_bl[block] = d3k3
                # Calculate integrated times
                int_tau_bl[block] = Integrate(time,cab)
                int_dc2_bl[block] = Integrate(time,dcab)
                int_d2c2_bl[block] = Integrate(time,d2cab)
                int_d3c2_bl[block] = Integrate(time,d3cab)
                int_ea_bl[block] = (1.0/int_tau_bl[block])*int_dc2_bl[block]
                # Calculate Components
                tau2_comp1_bl[block], tau2_comp2_bl[block] = Ea_Breakdown(A1,-dA1,(1.0/k1),int_tau_bl[block], dk1*(1.0/k1))
                taul_comp1_bl[block], taul_comp2_bl[block] = Ea_Breakdown(A2,-dA2,(1.0/k2),int_tau_bl[block], dk2*(1.0/k2))
                taui_comp1_bl[block], taui_comp2_bl[block] = Ea_Breakdown(A3,-dA3,(1.0/k3),int_tau_bl[block], dk3*(1.0/k3))
                comp_sum_bl[block] = tau2_comp1_bl[block] + tau2_comp2_bl[block] + taul_comp1_bl[block] + taul_comp2_bl[block] +taui_comp1_bl[block] + taui_comp2_bl[block]
                

                print("block <tau2> %s %s" % (int_tau_bl[block], int_ea_bl[block]))
                # Make predictions for the blocks
                for b in range(0,200):
                    beta[b]=0.5+float(b)*0.01
                    order2pred[b]=Order2_Prediction(beta[b],T,int_tau_bl[block], int_dc2_bl[block], int_d2c2_bl[block])
                    order3pred[b]=Order3_Prediction(beta[b],T,int_tau_bl[block], int_dc2_bl[block], int_d2c2_bl[block], int_d3c2_bl[block])
                order2pred_bl[block]=order2pred
                order3pred_bl[block]=order3pred
                
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
            err_d3A1 = Error(d3A1_bl)
            err_d3A2 = Error(d3A2_bl)
            err_d3A3 = Error(d3A3_bl)
            err_k1 = Error(k1_bl)
            err_k2 = Error(k2_bl)
            err_k3 = Error(k3_bl)
            err_dk1 = Error(dk1_bl)
            err_dk2 = Error(dk2_bl)
            err_dk3 = Error(dk3_bl)
            err_d2k1 = Error(d2k1_bl)
            err_d2k2 = Error(d2k2_bl)
            err_d2k3 = Error(d2k3_bl)
            err_d3k1 = Error(d3k1_bl)
            err_d3k2 = Error(d3k2_bl)
            err_d3k3 = Error(d3k3_bl)
            print dk1_bl
            print k1_bl
            err_int_tau = Error(int_tau_bl)
            err_int_dc2 = Error(int_dc2_bl)
            err_int_d2c2 = Error(int_d2c2_bl)
            err_int_d3c2 = Error(int_d3c2_bl)
            err_int_ea = Error(int_ea_bl)
            err_tau2_comp1 = Error(tau2_comp1_bl)
            err_tau2_comp2 = Error(tau2_comp2_bl)
            err_taul_comp1 = Error(taul_comp1_bl)
            err_taul_comp2 = Error(taul_comp2_bl)
            err_taui_comp1 = Error(taui_comp1_bl)
            err_taui_comp2 = Error(taui_comp2_bl)
            err_comp_sum = Error(comp_sum_bl)

            order2prederr=np.array(order2pred_bl).std(0)
            order2prederr=[x * t_val for x in order2prederr]
            order3prederr=np.array(order3pred_bl).std(0)
            order3prederr=[x * t_val for x in order3prederr]

            np.savetxt('fiterr_'+item1+'_'+str(mol_name)+'_c2.dat', np.c_[err_A1, err_A2, err_A3, err_k1, err_k2, err_k3, err_dA1, err_dA2, err_dA3, err_dk1, err_dk2, err_dk3], fmt='%s')
            np.savetxt('fiterr_'+item1+'_'+item2+'_'+str(mol_name)+'_c2.dat', np.c_[err_d2A1, err_d2A2, err_d2A3, err_d2k1, err_d2k2, err_d2k3], fmt='%s')
            np.savetxt('fiterr_'+item1+'_'+item2+'_'+item3+'_'+str(mol_name)+'_c2.dat', np.c_[err_d3A1, err_d3A2, err_d3A3, err_d3k1, err_d3k2, err_d3k3], fmt='%s')
            np.savetxt('fiterr_'+item1+'_'+item2+'_'+item3+'_'+str(mol_name)+'_intc2.dat', np.c_[err_int_tau, err_int_dc2, err_int_d2c2, err_int_d3c2, err_int_ea, err_tau2_comp1, err_tau2_comp2, err_taul_comp1, err_taul_comp2, err_taui_comp1, err_taui_comp2, err_comp_sum], fmt='%s')
            np.savetxt('fiterr_'+item1+'_'+item2+'_'+item3+'_'+str(mol_name)+'_errc2.dat', np.c_[order2prederr,order3prederr], fmt='%s')
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
            item3=item1
            # Read and Fit Corr
            time, cab, dcab = np.genfromtxt(item1 + '_' + str(mol_name) + '_c2.dat', usecols=(0,1,3), unpack=True)
            popt_cab, pcov_cab = curve_fit(Corr_Fit, time[:cut], cab[:cut], p0=(.2,.5,.5,.025,.45,3.1),bounds=((0, 0, 0, 0, 0, 0), (1, 1, 1, np.inf, np.inf, np.inf)))
            A=popt_cab[:3]
            k=popt_cab[3:]
            ksrt=sorted(k)
            A1=A[np.where(k==ksrt[0])]
            A2=A[np.where(k==ksrt[1])]
            A3=1-A1-A2
            k1=ksrt[0]
            k2=ksrt[1]
            k3=ksrt[2]
            # Fit Alt Correlation 
            popt_altcab, pcov_altcab = curve_fit(Alt_CorrFit, time[:cut], cab[:cut], p0=(.77799,.14109,.09,.483,2.5,40.,1.0),bounds=((0, 0, 0,0, .1, 0, 0), (1, 1, 1, np.inf,np.inf, np.inf,np.inf)))
            A_ALT=[popt_altcab[0],popt_altcab[1],popt_altcab[2]]
            k_ALT=[popt_altcab[3],popt_altcab[4],popt_altcab[5]]
            Y_ALT=popt_altcab[6]
            aA1=A_ALT[0]
            aA2=A_ALT[1]
            aA3=A_ALT[2]
            ak1=k_ALT[0]
            ak2=k_ALT[1]
            ak3=k_ALT[2]
            Y=Y_ALT
            # Fit Alt Derivative
            popt_altdcab, pcov_altdcab = curve_fit(Alt_dCorrFit, time[:cut], dcab[:cut], p0=(0.2,0.3,0.4,1.0,2.0,3.0))

            # Fit First Derivative
            print("First Total Derivative:")
            popt_dcab, pcov_dcab = curve_fit(dCorr_Fit,time, dcab,sigma=sig_arr, p0=(0.2,0.3,0.4,1.0,2.0,3.0))
            dA=popt_dcab[:3]
            dk=popt_dcab[3:]
            dA1=dA[0]
            dA2=dA[1]
            dA3=0.0-dA1-dA2
            dk1=dk[0]
            dk2=dk[1]
            dk3=dk[2]
            # Fit Second Derivative
            print("Second Total Derivative:")
            d2cab =np.genfromtxt(item1 + '_' + item2 + '_' + str(mol_name) + '_c2.dat', usecols=(5), unpack=True)
            popt_d2cab, pcov_d2cab = curve_fit(d2Corr_Fit, time, d2cab, p0=(0.2,0.3,0.4,1.2,0.3,0.4))
            d2A=popt_d2cab[:3]
            d2k=popt_d2cab[3:]
            d2A1=d2A[0]
            d2A2=d2A[1]
            d2A3=d2A[2]
            d2k1=d2k[0]
            d2k2=d2k[1]
            d2k3=d2k[2]
            # Read and Fit Third Derivative
            print("Third Total Derivative:")
            d3cab = np.genfromtxt(item1+'_'+item2+'_'+item3+'_'+mol_name+'_c2.dat', usecols=(1), unpack=True)
            popt_d3cab, pcov_d3cab = curve_fit(d3Corr_Fit, time, d3cab, p0=(0.2,0.3,0.4,0.2,0.3,0.4))
            d3A=popt_d3cab[:3]
            d3k=popt_d3cab[3:]
            d3A1=d3A[0]
            d3A2=d3A[1]
            d3A3=d3A[2]
            d3k1=d3k[0]
            d3k2=d3k[1]
            d3k3=d3k[2]
            # Read in Uncertainties
            err_A1, err_A2, err_A3, err_k1, err_k2, err_k3, err_dA1, err_dA2, err_dA3, err_dk1, err_dk2, err_dk3 = np.genfromtxt('fiterr_'+item1+'_'+str(mol_name)+'_c2.dat', usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True)
            err_d2A1, err_d2A2, err_d2A3, err_d2k1, err_d2k2, err_d2k3 = np.genfromtxt('fiterr_'+item1+'_'+item2+'_'+str(mol_name)+'_c2.dat', usecols=(0,1,2,3,4,5), unpack=True)
            err_d3A1, err_d3A2, err_d3A3, err_d3k1, err_d3k2, err_d3k3 = np.genfromtxt('fiterr_'+item1+'_'+item2+'_'+item3+'_'+str(mol_name)+'_c2.dat', usecols=(0,1,2,3,4,5), unpack=True)
            err_int_tau, err_int_dc2, err_int_d2c2, err_int_d3c2, err_int_ea, err_tau2_comp1, err_tau2_comp2, err_taul_comp1, err_taul_comp2, err_taui_comp1, err_taui_comp2, err_comp_sum= np.genfromtxt('fiterr_'+item1+'_'+item2+'_'+item3+'_'+str(mol_name)+'_intc2.dat', usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True)
            order2prederr = np.genfromtxt('fiterr_'+item1+'_'+item2+'_'+item3+'_'+str(mol_name)+'_errc2.dat', usecols=(0), unpack=True)
            # Calculate integrated times
            int_tau=Integrate(time,cab)
            int_dc2=Integrate(time,dcab)
            int_d2c2=Integrate(time,d2cab)
            int_d3c2=Integrate(time,d3cab)
            # Make predictions
            for b in range(0,200):
                beta[b]=0.5+float(b)*0.01
                order2pred[b]=Order2_Prediction(beta[b],T,int_tau, int_dc2, int_d2c2)
                order3pred[b]=Order3_Prediction(beta[b],T,int_tau, int_dc2, int_d2c2, int_d3c2)
            # Write out to log
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
            fout.write("d3A1: %s   err: %s\n" % (d3A1, err_d3A1))
            fout.write("d3A2: %s   err: %s\n" % (d3A2, err_d3A2))
            fout.write("d3A3: %s   err: %s\n" % (d3A3, err_d3A3))
            fout.write("d3k1: %s   err: %s\n" % (d3k1, err_d3k1))
            fout.write("d3k2: %s   err: %s\n" % (d3k2, err_d3k2))
            fout.write("d3k3: %s   err: %s\n" % (d3k3, err_d3k3))
            fout.write("\n")
            fout.write(" Ea Tau 2: %s   err: %s\n" % (dk1*(1.0/k1), E_Divi_Prop(dk1,err_dk1, k1, err_k1)))
            fout.write(" Tau 2 Time: %s +/- %s ps\n" % ((1.0/k1),(err_k1/k1*(1.0/k1))))
            fout.write(" Ea Librational: %s   err: %s\n" % (dk2*(1.0/k2), E_Divi_Prop(dk2,err_dk2, k2, err_k2)))
            fout.write(" Librational Time: %s +/- %s ps\n" % ((1.0/k2),(err_k2/k2*(1.0/k2)))) 
            fout.write(" Ea Inertial: %s   err: %s\n" % (dk3*(1.0/k3), E_Divi_Prop(dk3,err_dk3, k3, err_k3)))
            fout.write(" Inertial Time: %s +/- %s ps\n" % ((1.0/k3), (err_k3/k3*(1.0/k3))))
            fout.write("Integrated Times\n")
            fout.write("<tau2> = %s +/- %s\n" % (int_tau,err_int_tau))
            fout.write("intdc2 = %s +/- %s\n" % (int_dc2,err_int_dc2))
            fout.write("intd2c2 = %s +/- %s\n" % (int_d2c2,err_int_d2c2))
            fout.write("intd3c2 = %s +/- %s\n" % (int_d3c2, err_int_d3c2))
            fout.write("<EaTau2> = %s +/- %s\n" % (int_dc2*(-1.0/int_tau), err_int_ea))
            fout.write("---split Ea terms---\n")
            tau2_comp1, tau2_comp2 = Ea_Breakdown(A1,-dA1,(1.0/k1),int_tau, dk1*(1.0/k1))
            taul_comp1, taul_comp2 = Ea_Breakdown(A2,-dA2,(1.0/k2),int_tau, dk2*(1.0/k2))
            taui_comp1, taui_comp2 = Ea_Breakdown(A3,-dA3,(1.0/k3),int_tau, dk3*(1.0/k3))
            fout.write("tau2: %s +/- %s\n" % (tau2_comp1, err_tau2_comp1))
            fout.write("tau2: %s +/- %s\n" % (tau2_comp2[0], err_tau2_comp2))
            fout.write("taulib: %s +/- %s\n" % (taul_comp1, err_taul_comp1))
            fout.write("taulib: %s +/- %s\n" % (taul_comp2[0], err_taul_comp2))
            fout.write("tauiner: %s +/- %s\n" % (taui_comp1, err_taui_comp1))
            fout.write("tauiner: %s +/- %s\n" % (taui_comp2[0], err_taui_comp2))
            fout.write("sum: %s +/- %s\n" %((tau2_comp1+tau2_comp2+taul_comp1+taul_comp2+taui_comp1+taui_comp2),err_comp_sum))
            fout.write("-------Alt Cab---------\n")
            fout.write("tau iner %s\n" % (1.0/k_ALT[2]))
            fout.write("tau lib  %s\n" % (1.0/k_ALT[1]))
            fout.write("tau 2    %s\n" % (1.0/k_ALT[0]))
            fout.write("Y    %s\n" % (Y_ALT))

            
            fout.write("End Correlation function\n")
            fout.close()
            cabfit=Corr_Fit(time,A1,A2,A3,k1,k2,k3)
            altcabfit=Alt_CorrFit(time,*popt_altcab)
            altdcabfit=Alt_dCorrFit(time,*popt_altdcab)
            dcabfit=dCorr_Fit(time, dA1,dA2,dA3,dk1,dk2,dk3)
            d2cabfit=d2Corr_Fit(time, *popt_d2cab)
            d3cabfit=d3Corr_Fit(time, *popt_d3cab)
            # Save files
            np.savetxt('total_fit_'+str(item1)+'_'+str(item2)+'_'+str(mol_name)+'_c2.dat', np.c_[time,cabfit, dcabfit,d2cabfit], fmt='%s')
            np.savetxt('total_altfit_'+str(item1)+'_'+str(item2)+'_'+str(mol_name)+'_c2.dat', np.c_[time,altcabfit,altdcabfit], fmt='%s')
            np.savetxt('total_fit_'+str(item1)+'_'+str(item2)+'_'+str(item3)+'_'+str(mol_name)+'_c2.dat', np.c_[time, d3cabfit], fmt='%s')
            np.savetxt('int_fit_'+str(item1)+'_'+str(item2)+'_'+str(mol_name)+'_c2.dat', np.c_[beta, order2pred, order2prederr], fmt='%s')
            np.savetxt('int_fit_'+item1+'_'+item2+'_'+item3+'_'+mol_name+'_c2.dat', np.c_[beta, order3pred, order3prederr], fmt='%s')
        jindex+=1
    iindex+=1
