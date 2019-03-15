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
import numpy as np
import sys
from read_input import input

machine = str(sys.argv[1])
lines=[]
with open("Direct_Calculation_Fluctuations/src/dependencies/"+machine+"_header.dat",'r') as mfile:
    for line in mfile:
        lines.append(line)

# Calls the read_input class
inputparam = input("input_file")

njobs = int(np.ceil(inputparam.num_files/float(inputparam.num_rpj)))
# Generates a single job array that runs 50 jobs each.
dcn_file = "run_array.sh"
dcn = open(dcn_file, 'w')
for item in lines:
    dcn.write(item)
dcn.write("#SBATCH --array 0-%s\n" % (njobs-1))
dcn.write("\n")
dcn.write("cd $SLURM_SUBMIT_DIR\n")
if machine == "CRC":
    if inputparam.prog == "LAMMPS":
        dcn.write("module load lammps/22Aug2018\n")
    elif inputparam.prog == "CP2K":
        dcn.write("module load cp2k/6.0/popt\n")
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
    dcn.write("    mpirun lmp_mpi < in.nve -screen none\n")
elif inputparam.prog == "CP2K":
    if inputparam.cab == "IONPAIRING":
        dcn.write("    python ../../vel_reselect.py 6.94 18.998 298.15 restart.$CUR\n")
        dcn.write("    sed -i -e 's@vel.file@vel_'$CUR'_%s.vxyz@g' in.nve.cp2k\n" % inputparam.molec[0])
    dcn.write("    sed -i -e 's@restart.file@RESTART/restart.'$CUR'@g' in.nve.cp2k\n")
    dcn.write("    sed -i -e 's@traj.file0@traj_'$CUR'_%s.xyz@g' in.nve.cp2k\n" % inputparam.molec[0])
    if machine == "CRC":
        dcn.write("    mpirun cp2k.popt in.nve.cp2k\n")
    elif machine == "CORI":
        dcn.write("srun -n 128 -c 8 --cpu_bind=cores  /global/homes/c/cjmundy/SOURCEFORGE_CP2K/cp2k-5.1/cp2k/exe/CRAY-cori-Haswell-gnu/cp2k.psmp in.nve.cp2k > cp2k.out\n")

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
    dcn.write('    echo %s > mol.info\n' % inputparam.molec[0])
    dcn.write('    python ../../set_msd_calcs.py \n')
    dcn.write('    ../../flux_side\n\n')
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
if machine == "CRC":
    if inputparam.prog == "LAMMPS":
        nve.write('module load lammps/22Aug2018\n\n')
    elif inputparam.prog == "CP2K":
        nve.write('module load cp2k/6.0/popt\n\n')
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
    nve.write('mpirun lmp_mpi < in.nve -screen none\n\n\n')
elif inputparam.prog == "CP2K":
    if inputparam.cab == "IONPAIRING":
        nve.write("python ../../vel_reselect.py 6.94 18.998 298.15 restart.$CUR\n")
    nve.write("sed -i -e 's@restart.file@../../RESTART/restart.AAA@g' in.nve.cp2k\n")
    nve.write("sed -i -e 's@traj.file0@traj_AAA_%s.xyz@g' in.nve.cp2k\n" % inputparam.molec[0])
    if machine == "CRC":
        nve.write("    mpirun cp2k.popt in.nve.cp2k\n")
    elif machine == "CORI":
        nve.write("srun -n 128 -c 8 --cpu_bind=cores  /global/homes/c/cjmundy/SOURCEFORGE_CP2K/cp2k-5.1/cp2k/exe/CRAY-cori-Haswell-gnu/cp2k.psmp in.nve.cp2k > cp2k.out\n")
                                        
if inputparam.cab == "TRANSPORT":
    for i in range(0,inputparam.num_molecs):
        nve.write('echo %s > mol.info\n' % inputparam.molec[i])
        nve.write('python ../../set_msd_calcs.py\n')
        nve.write('../../msd_rot_calc < corr_calc.in\n\n')
    if inputparam.prog == "LAMMPS":
        nve.write('python ../../grab_press.py\n\n')
        nve.write('../../visc_calc\n')
elif inputparam.cab == "IONPAIRING":
    nve.write('echo %s > mol.info\n' % inputparam.molec[0])
    nve.write('python ../../set_msd_calcs.py \n')
    nve.write('../../flux_side\n\n')

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
iarr.write("#SBATCH --partition=sixhour\n")
iarr.write("#SBATCH --output=init_array%A_%a.out\n")
iarr.write("#SBATCH --nodes=1\n")
iarr.write("#SBATCH --ntasks-per-node=2\n")
iarr.write("#SBATCH --constraint=intel\n")
iarr.write("#SBATCH --mem=30G\n")
iarr.write("#SBATCH --time=06:00:00\n")
iarr.write("#SBATCH --array=0-%s\n" % (sep-1))


iarr.write("cd $SLURM_SUBMIT_DIR\n")

if inputparam.cab == "TRANSPORT":
    for i in range(inputparam.num_molecs):
        iarr.write("python init_segments.py $SLURM_ARRAY_TASK_ID flucts.inp msd %s\n" % inputparam.molec[i])
        iarr.write("python init_segments.py $SLURM_ARRAY_TASK_ID flucts.inp c2 %s\n" % inputparam.molec[i])
elif inputparam.cab == "IONPAIRING":
    iarr.write("python init_segments.py $SLURM_ARRAY_TASK_ID flucts.inp fsc_f water\n")
    iarr.write("python init_segments.py $SLURM_ARRAY_TASK_ID flucts.inp fsc_b water\n")

if inputparam.prog == "LAMMPS":
    iarr.write("python init_segments.py $SLURM_ARRAY_TASK_ID flucts.inp shear water\n")

iarr.close()

# Gen Combine Segments Array
# This submits the combine segmentspython code which does the final analysis.

do_file="combine_segments.sh"
darr = open(do_file, 'w')
darr.write("#!/bin/bash\n")
darr.write("#SBATCH --job-name=do_flucts\n")
darr.write("#SBATCH --partition=sixhour\n")
darr.write("#SBATCH --output=do_flucts.out\n")
darr.write("#SBATCH --nodes=1\n")
darr.write("#SBATCH --ntasks-per-node=20\n")
darr.write("#SBATCH --constraint=intel\n")
darr.write("#SBATCH --mem=100G\n")
darr.write("#SBATCH --time=06:00:00\n")
if inputparam.cab == "TRANSPORT":
    for i in range(inputparam.num_molecs):
        darr.write("python combine_segments.py flucts.inp msd %s\n" % (inputparam.molec[i]))
        darr.write("python combine_segments.py flucts.inp c2 %s\n" % (inputparam.molec[i]))
elif inputparam.cab == "IONPAIRING":
    darr.write("python combine_segments.py flucts.inp fsc_f %s\n" % (inputparam.molec[0]))
    darr.write("python combine_segments.py flucts.inp fsc_b %s\n" % (inputparam.molec[0]))

if inputparam.prog == "LAMMPS":
    darr.write("python combine_segments.py flucts.inp shear %s\n" % (inputparam.molec[0]))
darr.close()
