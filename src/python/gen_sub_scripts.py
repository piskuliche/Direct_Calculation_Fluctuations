# This python script generates the needed sub.sh, nve.sh, and job_array.sh scripts in the simulation directory i.e. /sim/
# required input: input_file
# produced output: setup_array.sh, nve.sh, job_array.sh
import numpy as np
from read_input import input

# Calls the read_input class
inputparam = input("input_file")

narrays = int(np.ceil(inputparam.num_files/5000.))
# Generates setupfiles submission array
# This section is the one that creates the job array that builds the filesystem - should be quick!
for i in range(narrays):
    start = int(i*5000)
    end = int(start + 5000.)
    sarr_file="setup_array"+str(i)+".sh"
    sarr=open(sarr_file, 'w')
    sarr.write('#MSUB -N setup_file_array\n')
    sarr.write('#MSUB -q sixhour\n')
    sarr.write('#MSUB -j oe\n')
    sarr.write('#MSUB -d ./\n')
    sarr.write('#MSUB -l nodes=1:ppn=1:intel,mem=5gb,walltime=6:00:00\n')
    sarr.write('#MSUB -t %s-%s\n\n\n' % (start,end))

    sarr.write('SEP=%s\n' % (inputparam.sep_config))
    sarr.write('START=%s\n' % (inputparam.start_config))
    sarr.write('CUR=$(( MOAB_JOBARRAYINDEX*SEP - SEP + START ))\n')
    sarr.write('cd $PBS_O_WORKDIR\n')

    sarr.write('mkdir FILES/$CUR\n')
    sarr.write('cp in.nve FILES/$CUR\n')
    sarr.write('cp nve.sh FILES/$CUR\n')
    sarr.write('cp msd_rot_calc FILES/$CUR\n')
    sarr.write('cp flux_side FILES/$CUR\n')
    sarr.write('cp grab_press.py FILES/$CUR\n')
    sarr.write('cp visc_calc FILES/$CUR\n')
    sarr.write('cp time.dat FILES/$CUR\n')
    sarr.write('cp set_msd_calcs.py FILES/$CUR\n')
    sarr.write('cp read_input.py FILES/$CUR\n')
    sarr.write('cp RESTART/restart.$CUR FILES/$CUR\n')
    sarr.write('cd FILES/$CUR\n')
    if inputparam.prog == "LAMMPS":
        sarr.write("sed -i -e 's@direct_calc_nve@nve_\'$CUR\'@g' nve.sh\n")
        sarr.write("sed -i -e 's@restart.file@../../RESTART/restart.\'$CUR\'@g' in.nve\n")
        for j in range(0, inputparam.num_molecs):
            sarr.write("sed -i -e 's@traj.file%s@traj_\'$CUR\'_%s.xyz@g' in.nve\n" % (str(j), inputparam.molec[j]))
            if inputparam.cab == 'IONPAIRING':
                sarr.write("sed -i -e 's@vel.file@vel_\'$CUR\'_%s.vxyz@g' in.nve\n" % inputparam.molec[j])
    elif inputparam.prog == "CP2K":
        sarr.write("sed -i -e 's@direct_calc_nve@nve_\'$CUR\'@g' nve.sh\n")
        sarr.write("sed -i -e 's@restart.file@../../RESTART/restart.\'$CUR\'@g' in.nve.cp2k\n")
        for j in range(0, inputparam.num_molecs):
            sarr.write("sed -i -e 's@traj.file%s@traj_\'$CUR\'_%s.xyz@g' in.nve.cp2k\n" % (str(j), inputparam.molec[j]))
            if inputparam.cab == 'IONPAIRING':
                sarr.write("sed -i -e 's@vel.file@vel_\'$CUR\'_%s.vxyz@g' in.nve.cp2k\n" % inputparam.molec[j])
    else: 
        print("Error: Incorrect input program in input_file")
    sarr.write('cd ../../\n')
    sarr.close()

# generate submit file
submitfile="setup_files"
sfi=open(submitfile,'w')
for i in range(narrays):
    sfi.write('msub setup_array%s.sh\n' % i)
    if i != narrays-1:
        sfi.write('sleep 20m\n')
    else:
        sfi.write('touch ../.flag_setcomplete')
sfi.close()

# Generate Job Array File
ja_file="job_array.sh"
ja=open(ja_file, 'w')

ja.write('#MSUB -N direct_calc_nve\n')
ja.write('#MSUB -q sixhour\n')
ja.write('#MSUB -j oe\n')
ja.write('#MSUB -d ./\n')
ja.write('#MSUB -l nodes=1:ppn=2:intel,mem=5gb,walltime=6:00:00\n')
ja.write('#MSUB -t AAA-BBB\n\n\n')

ja.write('SEP=%s\n' % (inputparam.sep_config))
ja.write('START=%s\n' % (inputparam.start_config))
ja.write('CUR=$(( MOAB_JOBARRAYINDEX*SEP - SEP + START ))\n')
ja.write('cd $PBS_O_WORKDIR\n')
ja.write('cd FILES/$CUR\n\n\n')

if inputparam.prog == "LAMMPS":
    ja.write('module load lammps/11Aug17\n\n')
elif inputparam.prog == "CP2K":
    ja.write('module load cp2k/6.0/popt\n\n')

ja.write('echo Time is `date` > array_$MOAB_JOBARRAYINDEX.o\n')
ja.write('echo Directory is `pwd` >> array_$MOAB_JOBARRAYINDEX.o\n\n\n')
if inputparam.prog == "LAMMPS":
    ja.write('mpirun lmp_mpi < in.nve -screen none\n\n\n')
elif inputparam.prog == "CP2K":
    ja.write('mpirun -np 2 cp2k.popt in.nve.cp2k \n\n\n')
if inputparam.cab == "TRANSPORT":
    for i in range(0,inputparam.num_molecs):
        ja.write('echo %s > mol.info\n' % inputparam.molec[i])
        ja.write('python ../../set_msd_calcs.py \n')
        ja.write('./msd_rot_calc < msd_rot_calc.in\n\n')

    if inputparam.prog == "LAMMPS":
        ja.write('python grab_press.py\n')
        ja.write('./visc_calc\n\n')
elif inputparam.cab == "IONPAIRING":
    ja.write('echo %s > mol.info\n' % inputparam.molec[0])
    ja.write('python ../../set_msd_calcs.py \n')
    ja.write('./flux_side\n\n')


ja.write('echo Ending Time is `date` >> array_$MOAB_JOBARRAYINDEX.o\n')

# Generate NVE Input File
nve_file="nve.sh"
nve=open(nve_file, 'w')

nve.write('#MSUB -N direct_calc_nve\n')
nve.write('#MSUB -q sixhour\n')
nve.write('#MSUB -j oe\n')
nve.write('#MSUB -d ./\n')
nve.write('#MSUB -l nodes=1:ppn=10:intel,mem=5gb,walltime=6:00:00\n\n\n')


if inputparam.prog == "LAMMPS":
    nve.write('module load lammps/11Aug17\n\n')
elif inputparam.prog == "CP2K":
    nve.write('module load cp2k/6.0/popt\n\n')
            
nve.write('echo Time is `date` > array_$MOAB_JOBARRAYINDEX.o\n')
nve.write('echo Directory is `pwd` >> array_$MOAB_JOBARRAYINDEX.o\n\n\n')

if inputparam.prog == "LAMMPS":
    nve.write('mpirun lmp_mpi < in.nve -screen none\n\n\n')
elif inputparam.prog == "CP2K":
    nve.write('mpirun -np 2 cp2k.popt in.nve.cp2k \n\n\n')
if inputparam.cab == "TRANSPORT":
    for i in range(0,inputparam.num_molecs):
        nve.write('echo %s > mol.info\n' % inputparam.molec[i])
        nve.write('python ../../set_msd_calcs.py\n')
        nve.write('./msd_rot_calc < msd_rot_calc.in\n\n')
    if inputparam.prog == "LAMMPS":
        nve.write('python grab_press.py\n\n')
        nve.write('./visc_calc\n')
elif inputparam.cab == "IONPAIRING":
    nve.write('echo %s > mol.info\n' % inputparam.molec[0])
    nve.write('python ../../set_msd_calcs.py \n')
    nve.write('./flux_side\n\n')

nve.write('echo Ending Time is `date` >> array_$MOAB_JOBARRAYINDEX.o\n')


nve.close()
