# This python script generates the needed sub.sh, nve.sh, and job_array.sh scripts in the simulation directory i.e. /sim/
# required input: input_file
# produced output: sub.sh, nve.sh, job_array.sh

from read_input import input

# Calls the read_input class
inputparam = input("input_file")


# Generates setupfiles submission script

sub_file="sub.sh"
sub=open(sub_file, 'w')

sub.write('#MSUB -N setupfiles\n')
sub.write('#MSUB -q sixhour\n')
sub.write('#MSUB -d ./\n')
sub.write('#MSUB -j oe\n')
sub.write('#MSUB -l nodes=1:ppn=10:intel,mem=5gb,walltime=6:00:00\n\n')

sub.write('bash setup_files\n\n')

sub.close()

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
