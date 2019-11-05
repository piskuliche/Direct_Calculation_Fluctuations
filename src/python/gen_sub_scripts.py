# This python script generates the needed sub.sh, nve.sh, and job_array.sh scripts in the simulation directory i.e. /sim/
"""
This python script generates the needed files in the simulation directory.
Requires an input_file
This produces the following output:
    run_array.sh
    nve.sh
    init_segments.sh
    combine_segments.sh
"""

# These set global variables and should be changed based on your system
lammpsload="module load lammps/12Dec2018"
lammpsruncmd="mpirun lmp_mpi < in.nve -screen none"
cp2kload="module load cp2k/6.0/popt"
cp2kruncmd="    mpirun cp2k.popt in.nve.cp2k"
partition="sixhour"
#cp2kruncmd="    srun -n 128 -c 8 --cpu_bind=cores  /global/homes/c/cjmundy/SOURCEFORGE_CP2K/cp2k-5.1/cp2k/exe/CRAY-cori-Haswell-gnu/cp2k.psmp in.nve.cp2k > cp2k.out"


def read_machine():
    f = open('machine.name','r')
    machine = str(f.readline().strip())
    path = str(f.readline().strip())
    return machine,path


import numpy as np
import sys
from src.python.read_input import user_input

machine, homepath = read_machine()

lines=[]
with open(homepath+"/src/dependencies/"+machine+"_header.dat",'r') as mfile:
    for line in mfile:
        lines.append(line)

# Calls the read_input class
inputparam = user_input("input_file")

njobs = int(np.ceil(inputparam.num_files/float(inputparam.num_rpj)))
# Generates a single job array that runs 50 jobs each.
dcn_file = "run_array.sh"
dcn = open(dcn_file, 'w')
for item in lines:
    dcn.write(item)
dcn.write("#SBATCH --array 0-%s\n" % (njobs-1))
dcn.write("\n")
dcn.write("module load Dir_Calc_Fluct\n")
dcn.write("cd $SLURM_SUBMIT_DIR\n")
if inputparam.prog == "LAMMPS":
    dcn.write("%s\n" % lammpsload)
elif inputparam.prog == "CP2K":
    dcn.write("%s\n" % cp2kload)
else:
    dcn.write("Error: Incorrect program in input_file\n")
    exit(1)
elif machine == "CORI":
    dcn.write("export OMP_NUM_THREADS=8\n")
    dcn.write("export OMP_PLACES=threads\n")
    dcn.write("export OMP_PROC_BIND=spread\n")

dcn.write("\n")
dcn.write("SEP=%s\n" % int(inputparam.sep_config))
dcn.write("START=%s\n" % int(inputparam.start_config))
dcn.write("\n")
dcn.write("START_JOBS=1\n")
dcn.write("NUM_RPJ=%s\n" % inputparam.num_rpj)
dcn.write("\n")
dcn.write("for (( N = $START_JOBS; N <= $NUM_RPJ; N++ ))\n")
dcn.write("do\n")
dcn.write("    CUR=$(( START + N*SEP + SEP*NUM_RPJ*SLURM_ARRAY_TASK_ID - SEP))\n")
dcn.write("    echo $CUR\n")
dcn.write("    mkdir FILES/$CUR\n")
if inputparam.prog == "LAMMPS":
    dcn.write("    cp in.nve nve.sh FILES/$CUR/\n")
elif inputparam.prog == "CP2K":
    dcn.write("    cp in.nve.cp2k nve.sh FILES/$CUR/\n")
else:
    dcn.write("Error: Incorrect program in input_file")
dcn.write("    cp RESTART/restart.$CUR FILES/$CUR\n")
dcn.write("    cd FILES/$CUR\n")
dcn.write("    echo Time is `date` > array_$SLURM_ARRAY_TASK_ID.o\n")
dcn.write("    echo Directory is `pwd` >> array_$SLURM_ARRAY_TASK_ID.o\n")
dcn.write("    \n\n")
dcn.write("    sed -i -e 's@direct_calc_nve@nve_'$CUR'@g' nve.sh\n")
if inputparam.prog == "LAMMPS":
    dcn.write("    sed -i -e 's@restart.file@../../RESTART/restart.'$CUR'@g' in.nve\n")
    for j in range(0,inputparam.num_molecs):
        dcn.write("    sed -i -e 's@traj.file%s@traj_'$CUR'_%s.xyz@g' in.nve\n" % (str(j),inputparam.molec[j]))
    dcn.write("%s\n" % lammpsruncmd)
elif inputparam.prog == "CP2K":
    if inputparam.cab == "IONPAIRING":
        dcn.write("    vel_reselect.py 6.94 18.998 298.15 restart.$CUR\n")
        dcn.write("    sed -i -e 's@vel.file@vel_'$CUR'_%s.vxyz@g' in.nve.cp2k\n" % inputparam.molec[0])
    dcn.write("    sed -i -e 's@restart.file@RESTART/restart.'$CUR'@g' in.nve.cp2k\n")
    dcn.write("    sed -i -e 's@traj.file0@traj_'$CUR'_%s.xyz@g' in.nve.cp2k\n" % inputparam.molec[0])
    dcn.write("%s\n" % runcmd)

dcn.write("    \n\n")
if inputparam.cab == "TRANSPORT":
    for i in range(0,inputparam.num_molecs):
        dcn.write("    echo %s > mol.info\n" % inputparam.molec[i])
        dcn.write("    set_msd_calcs.py\n")
        if inputparam.molec[i] == "water":
            dcn.write("    msd_rot_calc < corr_calc.in\n")
            dcn.write("    jump_rot.py\n")
        else:
            dcn.write("    matom_msd_rot_calc < corr_calc.in\n")
    dcn.write("    \n")
    if inputparam.prog == "LAMMPS":
        dcn.write("    grab_press.py\n")
        dcn.write("    visc_calc\n")
        dcn.write("    \n")
elif inputparam.cab == "IONPAIRING":
    dcn.write('    echo %s > mol.info\n' % inputparam.molec[0])
    dcn.write('    set_msd_calcs.py \n')
    dcn.write('    flux_side\n\n')
dcn.write('    touch corr_complete\n')
dcn.write("    echo Ending Time is `date` >> array_$SLURM_ARRAY_TASK_ID.o\n")
dcn.write("    rm mol.info\n")
dcn.write("    rm traj*.xyz\n")
if inputparam.cab == "IONPAIRING":
    dcn.write("    rm vel*.vxyz\n")
dcn.write("    cd ../../\n")
dcn.write("done\n")
dcn.close()


# Generate NVE Input File
#   This is the file that is generated in every directory that needs it if
#   something doesn't run correctly.
nve_file="nve.sh"
nve=open(nve_file, 'w')

for item in lines:
    nve.write(item)
nve.write("module load Dir_Calc_Fluct\n")
if inputparam.prog == "LAMMPS":
    nve.write('%s\n\n' % lammpsload)
elif inputparam.prog == "CP2K":
    nve.write('%s\n\n' cp2kload)
elif machine == "CORI":
    nve.write("export OMP_NUM_THREADS=8\n")
    nve.write("export OMP_PLACES=threads\n")
    nve.write("export OMP_PROC_BIND=spread\n")
nve.write('echo Time is `date` > array.o\n')
nve.write('echo Directory is `pwd` >> array.o\n\n\n')

if inputparam.prog == "LAMMPS":
    nve.write("sed -i -e 's@restart.file@../../RESTART/restart.AAA@g' in.nve\n")
    for j in range(0,inputparam.num_molecs):
        nve.write("sed -i -e 's@traj.file%s@traj_AAA_%s.xyz@g' in.nve\n" % (str(j),inputparam.molec[j]))
    nve.write('%s\n\n\n' % lammpsruncmd)
elif inputparam.prog == "CP2K":
    if inputparam.cab == "IONPAIRING":
        nve.write("vel_reselect.py 6.94 18.998 298.15 restart.$CUR\n")
    nve.write("sed -i -e 's@restart.file@../../RESTART/restart.AAA@g' in.nve.cp2k\n")
    nve.write("sed -i -e 's@traj.file0@traj_AAA_%s.xyz@g' in.nve.cp2k\n" % inputparam.molec[0])
    nve.write("%s\n" % cp2kruncmd)
                                        
if inputparam.cab == "TRANSPORT":
    for i in range(0,inputparam.num_molecs):
        nve.write('echo %s > mol.info\n' % inputparam.molec[i])
        nve.write('set_msd_calcs.py\n')
        if inputparam.molec[i] == "water":
            nve.write('msd_rot_calc < corr_calc.in\n\n')
            nve.write('jump_rot.py\n')
        else:
            nve.write('matom_msd_rot_calc < corr_calc.in\n\n')
    if inputparam.prog == "LAMMPS":
        nve.write('grab_press.py\n\n')
        nve.write('visc_calc\n')
elif inputparam.cab == "IONPAIRING":
    nve.write('echo %s > mol.info\n' % inputparam.molec[0])
    nve.write('set_msd_calcs.py \n')
    nve.write('flux_side\n\n')

nve.write('touch corr_complete\n')
nve.write('echo Ending Time is `date` >> array.o\n')


nve.close()

# Gen Init Array
init_file="init_segments.sh"
iarr = open(init_file,'w')
sep = int(inputparam.num_files/inputparam.segsplit)
if sep % inputparam.nblocks != 0:
    print("Num files and Nblocks do not divide evenly, aborting")
    sys.exit(1)
iarr.write("#!/bin/bash\n")
iarr.write("#SBATCH --job-name=init_array\n")
iarr.write("#SBATCH --partition=%s\n" % partition)
iarr.write("#SBATCH --output=init_array%A_%a.out\n")
iarr.write("#SBATCH --nodes=1\n")
iarr.write("#SBATCH --ntasks-per-node=2\n")
iarr.write("#SBATCH --constraint=intel\n")
iarr.write("#SBATCH --mem=30G\n")
iarr.write("#SBATCH --time=06:00:00\n")
iarr.write("#SBATCH --array=0-%s\n" % (sep-1))

iarr.write("module load Dir_Calc_Fluct\n")

iarr.write("cd $SLURM_SUBMIT_DIR\n")

if inputparam.cab == "TRANSPORT":
    for i in range(inputparam.num_molecs):
        iarr.write("init_segments.py $SLURM_ARRAY_TASK_ID flucts.inp msd %s\n" % inputparam.molec[i])
        iarr.write("init_segments.py $SLURM_ARRAY_TASK_ID flucts.inp c2 %s\n" % inputparam.molec[i])
elif inputparam.cab == "IONPAIRING":
    iarr.write("init_segments.py $SLURM_ARRAY_TASK_ID flucts.inp fsc_f water\n")
    iarr.write("init_segments.py $SLURM_ARRAY_TASK_ID flucts.inp fsc_b water\n")

if inputparam.prog == "LAMMPS":
    iarr.write("init_segments.py $SLURM_ARRAY_TASK_ID flucts.inp shear water\n")

iarr.close()

# Gen Combine Segments Array
# This submits the combine segmentspython code which does the final analysis.

do_file="combine_segments.sh"
darr = open(do_file, 'w')
darr.write("#!/bin/bash\n")
darr.write("#SBATCH --job-name=do_flucts\n")
darr.write("#SBATCH --partition=%s\n" % partition)
darr.write("#SBATCH --output=do_flucts.out\n")
darr.write("#SBATCH --nodes=1\n")
darr.write("#SBATCH --ntasks-per-node=20\n")
darr.write("#SBATCH --constraint=intel\n")
darr.write("#SBATCH --mem=100G\n")
darr.write("#SBATCH --time=06:00:00\n")
darr.write("module load Dir_Calc_Fluct\n")
if inputparam.cab == "TRANSPORT":
    for i in range(inputparam.num_molecs):
        darr.write("combine_segments.py flucts.inp msd %s\n" % (inputparam.molec[i]))
        darr.write("combine_segments.py flucts.inp c2 %s\n" % (inputparam.molec[i]))
elif inputparam.cab == "IONPAIRING":
    darr.write("combine_segments.py flucts.inp fsc_f %s\n" % (inputparam.molec[0]))
    darr.write("combine_segments.py flucts.inp fsc_b %s\n" % (inputparam.molec[0]))

if inputparam.prog == "LAMMPS":
    darr.write("combine_segments.py flucts.inp shear %s\n" % (inputparam.molec[0]))
darr.close()
