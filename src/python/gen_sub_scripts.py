# This python script generates the needed sub.sh, nve.sh, and job_array.sh scripts in the simulation directory i.e. /sim/
# required input: input_file
# produced output: setup_array.sh, nve.sh, job_array.sh
"""
This python script generates the needed files in the simulation directory.
Requires an input_file
This produces the following output:
    run_array#.sh
    run_array
    nve.sh
    init_array.sh
    do_flucts
"""
import numpy as np
import sys
from read_input import input

# Calls the read_input class
inputparam = input("input_file")

njobs = int(np.ceil(inputparam.num_files/50.))
# Generates a single job array that runs 50 jobs each.
dcn_file = "run_array.sh"
dcn = open(dcn_file, 'w')
dcn.write("#MSUB -N direct_calc_nve\n")
dcn.write("#MSUB -q sixhour\n")
dcn.write("#MSUB -j oe\n")
dcn.write("#MSUB -d ./\n")
dcn.write("#MSUB -l nodes=1:ppn=2:intel,mem=10gb,walltime=6:00:00\n")
dcn.write("#MSUB -t 0-%s\n" % (njobs-1))
dcn.write("\n")
dcn.write("cd $PBS_O_WORKDIR\n")
if inputparam.prog == "LAMMPS":
    dcn.write("module load lammps/11Aug17\n")
elif inputparam.prog == "CP2K":
    dcn.write("module load cp2k/6.0/popt\n")
else:
    dcn.write("Error: Incorrect program in input_file\n")
    exit(1)
dcn.write("\n")
dcn.write("SEP=1000\n")
dcn.write("START=1001000\n")
dcn.write("\n")
dcn.write("START_JOBS=1\n")
dcn.write("NUM_RPJ=50\n")
dcn.write("\n")
dcn.write("for (( N = $START_JOBS; N <= $NUM_RPJ; N++ ))\n")
dcn.write("do\n")
dcn.write("    CUR=$(( START + N*SEP + SEP*NUM_RPJ*MOAB_JOBARRAYINDEX - SEP))\n")
dcn.write("    echo $CUR\n")
dcn.write("    mkdir FILES/$CUR\n")
dcn.write("    cp in.nve nve.sh FILES/$CUR/\n")
dcn.write("    cp RESTART/restart.$CUR FILES/$CUR\n")
dcn.write("    cd FILES/$CUR\n")
dcn.write("    echo Time is `date` > array_$MOAB_JOBARRAYINDEX.o\n")
dcn.write("    echo Directory is `pwd` >> array_$MOAB_JOBARRAYINDEX.o\n")
dcn.write("    \n\n")
dcn.write("    sed -i -e 's@direct_calc_nve@nve_'$CUR'@g' nve.sh\n")
if inputparam.prog == "LAMMPS":
    dcn.write("    sed -i -e 's@restart.file@../../RESTART/restart.'$CUR'@g' in.nve\n")
    for j in range(0,inputparam.num_molecs):
        dcn.write("    sed -i -e 's@traj.file%s@traj_'$CUR'_%s.xyz@g' in.nve\n" % (str(j),inputparam.molec[j]))
    dcn.write("    mpirun lmp_mpi < in.nve -screen none\n")
elif inputparam.prog == "CP2K":
    dcn.write("    sed -i -e 's@restart.file@../../RESTART/restart.'$CUR'@g' in.nve.cp2k\n")
    dcn.write("    sed -i -e 's@traj.file0@traj_'$CUR'_%s.xyz@g' in.nve.cp2k\n" % inputparam.molec[0])
    dcn.write("    mpirun cp2k.popt in.nve.cp2k\n")
dcn.write("    \n\n")
if inputparam.cab == "TRANSPORT":
    for i in range(0,inputparam.num_molecs):
        dcn.write("    echo %s > mol.info\n" % inputparam.molec[i])
        dcn.write("    python ../../set_msd_calcs.py\n")
        dcn.write("    ../../msd_rot_calc < corr_calc.in\n")
    dcn.write("    \n")
    if inputparam.prog == "LAMMPS":
        dcn.write("    python ../../grab_press.py\n")
        dcn.write("    ../../visc_calc\n")
        dcn.write("    \n")
elif inputparam.cab == "IONPAIRING":
    dcn.write('echo %s > mol.info\n' % inputparam.molec[0])
    dcn.write('python ../../set_msd_calcs.py \n')
    dcn.write('../../flux_side\n\n')
dcn.write("    echo Ending Time is `date` >> array_$MOAB_JOBARRAYINDEX.o\n")
dcn.write("    rm mol.info\n")
dcn.write("    rm traj*.xyz\n")
dcn.write("    cd ../../\n")
dcn.write("done\n")
dcn.close()


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
    nve.write("sed -i -e 's@restart.file@../../RESTART/restart.AAA@g' in.nve\n")
    for j in range(0,inputparam.num_molecs):
        nve.write("sed -i -e 's@traj.file%s@traj_AAA_%s.xyz@g' in.nve.cp2k\n" % (str(j),inputparam.molec[j]))
    nve.write('mpirun lmp_mpi < in.nve -screen none\n\n\n')
elif inputparam.prog == "CP2K":
    nve.write("sed -i -e 's@restart.file@../../RESTART/restart.AAA@g' in.nve.cp2k\n")
    nve.write("sed -i -e 's@traj.file0@traj_AAA_%s.xyz@g' in.nve.cp2k\n" % inputparam.molec[0])
    nve.write('mpirun -np 2 cp2k.popt in.nve.cp2k \n\n\n')
if inputparam.cab == "TRANSPORT":
    for i in range(0,inputparam.num_molecs):
        nve.write('echo %s > mol.info\n' % inputparam.molec[i])
        nve.write('python ../../set_msd_calcs.py\n')
        nve.write('../../msd_rot_calc < msd_rot_calc.in\n\n')
    if inputparam.prog == "LAMMPS":
        nve.write('python ../../grab_press.py\n\n')
        nve.write('../../visc_calc\n')
elif inputparam.cab == "IONPAIRING":
    nve.write('echo %s > mol.info\n' % inputparam.molec[0])
    nve.write('python ../../set_msd_calcs.py \n')
    nve.write('../../flux_side\n\n')

nve.write('echo Ending Time is `date` >> array_$MOAB_JOBARRAYINDEX.o\n')


nve.close()

# Gen Init Array
init_file="init_array.sh"
iarr = open(init_file,'w')
sep = int(inputparam.num_files/500.)
if sep % inputparam.nblocks != 0:
    print("Num files and Nblocks do not divide evenly, aborting")
    sys.exit(1)

iarr.write("#MSUB -N init_array\n")
iarr.write("#MSUB -q sixhour\n")
iarr.write("#MSUB -j oe\n")
iarr.write("#MSUB -d ./\n")
iarr.write("#MSUB -l nodes=1:ppn=2:intel,mem=30gb,walltime=6:00:00\n")
iarr.write("#MSUB -t 0-%s\n" % (sep-1))


iarr.write("cd $PBS_O_WORKDIR\n")

if inputparam.cab == "TRANSPORT":
    for i in range(inputparam.num_molecs):
        iarr.write("python init_flucts.py $MOAB_JOBARRAYINDEX flucts.inp msd %s\n" % inputparam.molec[i])
        iarr.write("python init_flucts.py $MOAB_JOBARRAYINDEX flucts.inp c2 %s\n" % inputparam.molec[i])
elif inputparam.cab == "IONPAIRING":
    iarr.write("python init_flucts.py $MOAB_JOBARRAYINDEX flucts.inp fsc water\n")

if inputparam.prog == "LAMMPS":
    iarr.write("python init_flucts.py $MOAB_JOBARRAYINDEX flucts.inp shear water\n")

iarr.close()

# Gen Do_Flucts Array

do_file="do_fluctsub.sh"
darr = open(do_file, 'w')
darr.write("#MSUB -N do_flucts\n")
darr.write("#MSUB -q sixhour\n")
darr.write("#MSUB -j oe\n")
darr.write("#MSUB -m ae\n")
darr.write("#MSUB -d ./\n")
darr.write("#MSUB -l nodes=1:ppn=20:intel,mem=100gb,walltime=6:00:00\n")
if inputparam.cab == "TRANSPORT":
    for i in range(inputparam.num_molecs):
        darr.write("python do_flucts.py flucts.inp msd %s %s\n" % (inputparam.molec[i],inputparam.nblocks))
        darr.write("python do_flucts.py flucts.inp c2 %s %s\n" % (inputparam.molec[i],inputparam.nblocks))
elif inputparam.cab == "IONPAIRING":
    darr.write("python do_flucts.py flucts.inp fsc %s %s\n" % (inputparam.molec[0],inputparam.nblocks))

if inputparam.prog == "LAMMPS":
    darr.write("python do_flucts.py flucts.inp shear %s %s\n" % (inputparam.molec[0],inputparam.nblocks))
darr.close()
