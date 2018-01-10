import numpy as np

files = np.genfromtxt("file_names", unpack=True)
flucts = np.genfromtxt("vol_init.out", unpack=True)

for i in range(0, len(files)):
    filepath=str("FILES/"+str(int(files[i]))+"/msd_rot_calc.in")
    print("Created %s" % (filepath))
    f=open(filepath, 'w')
    f.write("# Numeric File\n")
    f.write("%s\n" % (str(int(files[i]))))
    f.write("# Number_of_Times Sep_of_Times:\n")
    f.write("400 0.050\n")
    f.write("# Volume\n")
    f.write("%s\n" % (flucts[i]))
    f.close()
    
print "MSD Input Setup"
