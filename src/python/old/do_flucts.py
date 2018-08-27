import numpy as np
import sys
from scipy import stats
from read_input import input

def block_energy(energy, start, end, itemindex):
    tmp = 0
    for i in range(start*sep,end*sep):
        tmp = tmp + energy[itemindex][i]
    tmp = tmp/(end*sep-start*sep)
    return tmp

def block_array_err(barray, noblocks):
    sigma=[]
    for i in range(len(barray)):
        sigma.append(np.std(barray[i])*t_val)
    return sigma

#Import Input File Parameters
inputparam=input('input_file')


val=1
#Read in command line arguments
fname = str(sys.argv[1])
corr_name = str(sys.argv[2])
mol_name = str(sys.argv[3])


#Student T Value
t_val=stats.t.ppf(0.975,inputparam.nblocks-1)/np.sqrt(inputparam.nblocks)


#Calculate the Number of Segments
sep=500
num_segs = int(inputparam.num_files/float(sep))
segs_per_block = np.floor(num_segs/float(inputparam.nblocks))
print segs_per_block

inp_n, inp_c=np.genfromtxt(fname, usecols=(0,1), dtype=(str,int),unpack=True)
time = np.zeros(inputparam.num_times)
corr = np.zeros((len(inp_n),len(inp_n),num_segs, inputparam.num_times))
d1corr = np.zeros((len(inp_n),len(inp_n),num_segs, inputparam.num_times))
d2corr = np.zeros((len(inp_n),len(inp_n),num_segs, inputparam.num_times))
d3corr = np.zeros((len(inp_n),len(inp_n),num_segs, inputparam.num_times))
d4corr = np.zeros((len(inp_n),len(inp_n),num_segs, inputparam.num_times))

energy=[]


if val == '-h':
    print("This function needs a single integer input to run.")
    print("     Usage: python do_average.py fname corr_func mol_name")
    exit()
else:            
    item1count=0
    # Create Arrays of Results
    print("---Start Read---")
    for item1 in inp_n:
        print("The current Item is: %s" % item1)
        item2count=0
        for item2 in inp_n:
            for seg in range(num_segs):
                # 4D Arrays of shape [item1][item2][seg][time]
                time,corr[item1count][item2count][seg],d1corr[item1count][item2count][seg],d2corr[item1count][item2count][seg],d3corr[item1count][item2count][seg],d4corr[item1count][item2count][seg]=np.genfromtxt('SEG/seg_'+str(seg)+'_'+item1+'_'+item2+'_'+mol_name+'_'+corr_name+'.dat', usecols=(0,1,2,3,4,5), unpack=True)

            item2count+=1
        item1count+=1
    # Do blocking
    print("---End Read---")
    for item1 in inp_n:
        energy.append(np.genfromtxt(item1+'_init.out'))
    item1count=0
    for item1 in inp_n:
        print item1
        item2count=0
        for item2 in inp_n:
            item3=item1
            item4=item1
            # Zero arrays
            blockcorr=np.zeros((inputparam.nblocks,inputparam.num_times))
            blockd1corr=np.zeros((inputparam.nblocks,inputparam.num_times))
            blockd2corr=np.zeros((inputparam.nblocks,inputparam.num_times))
            blockd3corr=np.zeros((inputparam.nblocks,inputparam.num_times))
            blockd4corr=np.zeros((inputparam.nblocks,inputparam.num_times))
            blockea=np.zeros((inputparam.nblocks,inputparam.num_times))
            totcorr=np.zeros(inputparam.num_times)
            totd1corr=np.zeros(inputparam.num_times)
            totd2corr=np.zeros(inputparam.num_times)
            totd3corr=np.zeros(inputparam.num_times)
            totd4corr=np.zeros(inputparam.num_times)
            totea=np.zeros(inputparam.num_times)
            errcorr=np.zeros(inputparam.num_times)
            errd1corr=np.zeros(inputparam.num_times)
            errd2corr=np.zeros(inputparam.num_times)
            errd3corr=np.zeros(inputparam.num_times)
            errd4corr=np.zeros(inputparam.num_times)
            errea=np.zeros(inputparam.num_times)
            for block in range(inputparam.nblocks):
                bstart=int(block*segs_per_block)
                bend=int((block+1)*segs_per_block)
                bdist=bend-bstart
                e1av=block_energy(energy,bstart,bend,item1count)
                e2av=block_energy(energy,bstart,bend,item2count)
                e3av=e1av
                e4av=e1av
                for i in range(inputparam.num_times):
                    for seg in range(bstart,bend):
                        blockcorr[block][i]+=corr[item1count][item2count][seg][i]
                        blockd1corr[block][i]+=d1corr[item1count][item2count][seg][i]
                        blockd2corr[block][i]+=d2corr[item1count][item2count][seg][i]
                        blockd3corr[block][i]+=d3corr[item1count][item2count][seg][i]
                        blockd4corr[block][i]+=d4corr[item1count][item2count][seg][i]
                    blockcorr[block][i]=blockcorr[block][i]/float(bdist)
                    blockd1corr[block][i]=blockd1corr[block][i]/float(bdist)-e1av*blockcorr[block][i]
                    blockd2corr[block][i]=blockd2corr[block][i]/float(bdist)
                    blockd3corr[block][i]=blockd3corr[block][i]/float(bdist)
                    blockd4corr[block][i]=blockd4corr[block][i]/float(bdist)
                blockea[block]=blockd1corr[block]/blockcorr[block]
                np.savetxt("bl_"+str(block)+"_"+item1+"_"+mol_name+"_"+corr_name+".dat", np.c_[time, blockcorr[block], blockd1corr[block],blockea[block]], fmt='%s')
                np.savetxt("bl_"+str(block)+"_"+item1+"_"+item2+"_"+mol_name+"_"+corr_name+".dat", np.c_[time,blockd2corr[block],blockd2corr[block],blockd2corr[block]], fmt='%s')
                np.savetxt("bl_"+str(block)+"_"+item1+"_"+item2+"_"+item3+"_"+mol_name+'_'+corr_name+".dat", np.c_[time, blockd3corr[block]],fmt='%s')
                np.savetxt("bl_"+str(block)+"_"+item1+"_"+item2+"_"+item3+"_"+item4+"_"+mol_name+'_'+corr_name+".dat", np.c_[time, blockd4corr[block]],fmt='%s')
            errcorr=np.array(blockcorr).std(0)
            errcorr=[x * t_val for x in errcorr]
            errd1corr=np.array(blockd1corr).std(0)
            errd1corr=[x * t_val for x in errd1corr]
            errd2corr=np.array(blockd2corr).std(0)
            errd2corr=[x * t_val for x in errd2corr]
            errd3corr=np.array(blockd3corr).std(0)
            errd3corr=[x * t_val for x in errd3corr]
            errd4corr=np.array(blockd4corr).std(0)
            errd4corr=[x * t_val for x in errd4corr]
            errea=np.array(blockea).std(0)
            errea=[x * t_val for x in errea]
            seg_start=0
            e1av=block_energy(energy,seg_start,num_segs,item1count)
            print e1av
            print np.average(energy[item1count])
            e2av=block_energy(energy,seg_start,num_segs,item2count)
            for i in range(inputparam.num_times):
                for seg in range(num_segs):
                    totcorr[i]+=corr[item1count][item2count][seg][i]
                    totd1corr[i]+=d1corr[item1count][item2count][seg][i]
                    totd2corr[i]+=d2corr[item1count][item2count][seg][i]
                    totd3corr[i]+=d3corr[item1count][item2count][seg][i]
                    totd4corr[i]+=d4corr[item1count][item2count][seg][i]
                totcorr[i]=totcorr[i]/float(num_segs)
                totd1corr[i]=totd1corr[i]/float(num_segs)-e1av*totcorr[i]
                totd2corr[i]=totd2corr[i]/float(num_segs)
                totd3corr[i]=totd3corr[i]/float(num_segs)
                totd4corr[i]=totd4corr[i]/float(num_segs)
            totea=totd1corr/totcorr
            np.savetxt(item1+"_"+mol_name+"_"+corr_name+".dat", np.c_[time, totcorr, errcorr, totd1corr,errd1corr, totea, errea], fmt='%s')
            np.savetxt(item1+"_"+item2+"_"+mol_name+"_"+corr_name+".dat", np.c_[time, totd2corr, errd2corr, totd2corr, errd2corr, totd2corr, errd2corr], fmt='%s')
            np.savetxt(item1+"_"+item2+"_"+item3+"_"+mol_name+"_"+corr_name+".dat", np.c_[time, totd3corr, errd3corr], fmt='%s')
            np.savetxt(item1+"_"+item2+"_"+item3+"_"+item4+"_"+mol_name+"_"+corr_name+".dat", np.c_[time, totd4corr, errd4corr], fmt='%s')

            
            item2count+=1
        item1count+=1
    
