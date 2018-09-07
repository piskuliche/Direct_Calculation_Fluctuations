import numpy as np
import sys
from read_input import input


inputparam = input("input_file")
sep = 500

val = sys.argv[1]
fname = str(sys.argv[2])
corr_name = str(sys.argv[3])
mol_name = str(sys.argv[4])

# Read In Time
time = np.genfromtxt('real_time.dat', usecols=0)
time = [0] + time
time = time[:inputparam.num_times]

# Read In Fluctuations
inp_n, inp_c=np.genfromtxt(fname, usecols=(0,1), dtype=(str,int),unpack=True)

# Initialize Arrays
d1=[sep]
d2=[sep]
d3=[sep]
d4=[sep]
corr=[]
d1corr=[]
d2corr=[]
d3corr=[]
d4corr=[]

for i in range(inputparam.num_times):
    corr.append(0)
    d1corr.append(0)
    d2corr.append(0)
    d3corr.append(0)
    d4corr.append(0)

energy=[]
tcab=[]

print("The Argument Provided is %s:" % val)
if val == '-h':
    print("This function needs a single integer input to run.")
    print("     Usage: python init_average.py sep fname corr_func")
    exit()
else:
    itmp=int(val)
    print itmp
    inputparam = input("input_file")
    fstart = itmp*sep
    fend = (itmp+1)*sep
    print("Starting Trajectory: %s\nEnding Trajectory: %s" % (fstart,fend))
    # Read in Correlation Functions
    fnames = np.genfromtxt('file_names', dtype='string', unpack=True)
    print("Choosing %s of %s files" % (sep, len(fnames)))
    fcab=["FILES/"+str(s)+"/"+str(corr_name)+"_"+str(s)+"_"+str(mol_name)+".dat" for s in fnames[fstart:fend]]
    item1count=0
    for item1 in inp_n:
        energy.append(np.genfromtxt(item1+'_init.out'))
    for i in range(sep):
        tcab.append(np.genfromtxt(fcab[i],usecols=1,unpack=True))
    print len(tcab)
    for item1 in inp_n:
        print item1
        item2count=0
        for item2 in inp_n:
            corr=np.zeros((inputparam.num_times))
            d1corr=np.zeros((inputparam.num_times))
            d2corr=np.zeros((inputparam.num_times))
            d3corr=np.zeros((inputparam.num_times))
            d4corr=np.zeros((inputparam.num_times))
            # Calculate Average Energy
            e1av = 0.0
            e2av = 0.0
            for i in range(sep):
                enum  = fstart+i
                e1av   += energy[item1count][enum]
                e2av   += energy[item2count][enum]
            e1av /= float(sep)
            e2av /= float(sep)

            for i in range(sep):
                enum=fstart+i
                corr=corr+tcab[i]
                d1=energy[item1count][enum]-e1av
                d2=energy[item2count][enum]-e2av
                d1corr=d1corr+tcab[i]*d1
                d2corr=d2corr+tcab[i]*d1*d2
                d3corr=d3corr+tcab[i]*d1*d2*d2
                d4corr=d4corr+tcab[i]*d1*d2*d2*d2
            corr[:] = [x / float(sep) for x in corr]
            d1corr[:] = [x / float(sep) for x in d1corr]
            d2corr[:] = [x / float(sep) for x in d2corr]
            d3corr[:] = [x / float(sep) for x in d3corr]
            d4corr[:] = [x / float(sep) for x in d4corr]
            print len(corr), len(time), len(d1corr), len(d2corr), len(d3corr), len(d4corr)
            np.savetxt('SEG/seg_'+str(int(val))+'_'+item1+'_'+item2+'_'+mol_name+'_'+corr_name+'.dat', np.c_[time, corr, d1corr, d2corr, d3corr, d4corr])
            item2count+=1
        item1count+=1

