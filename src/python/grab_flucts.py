import numpy as np

filenames='file_names'
num_file=0
avg=[]
fnames=[]
ke_list= np.array([])
v_list=np.array([])
lj_list=np.array([])
coul_list=np.array([])
ewald_list=np.array([])
e_list=np.array([])
vol_list=np.array([])
with open(filenames) as f:
    for l in f:
        filename='FILES/'+l.rstrip()+'/log.lammps'
        print filename
        lookup = 'Step'
        lookup2= 'Loop time'
        totlines=0
        with open(filename) as myFile:
            for num, line in enumerate(myFile, 1):
                if lookup in line:
                    startskip=num
                if lookup2 in line:
                    endskip=num
                totlines=num
        endskip=totlines-endskip-4
        e, v, ke, lj, coul, ewald,vol = np.genfromtxt(filename, skip_header=startskip,skip_footer=endskip,usecols=(2,4,5,6,7,10,11), unpack=True)
        ke_list = np.append(ke_list,ke[0])
        v_list = np.append(v_list,v[0])
        e_list = np.append(e_list,e[0])
        lj_list = np.append(lj_list,lj[0])
        coul_list = np.append(coul_list,coul[0])
        ewald_list = np.append(ewald_list,ewald[0])
        vol_list = np.append(vol_list, vol[0])

np.savetxt('ke_init.out', ke_list, fmt=['%.4f'])
np.savetxt('v_init.out', v_list, fmt=['%.4f'])
np.savetxt('e_init.out', e_list, fmt=['%.4f'])
np.savetxt('lj_init.out', lj_list, fmt=['%.4f'])
np.savetxt('coul_init.out', coul_list, fmt=['%.4f'])
np.savetxt('ewald_init.out', ewald_list, fmt=['%.4f'])
np.savetxt('vol_init.out', vol_list, fmt=['%.4f'])
