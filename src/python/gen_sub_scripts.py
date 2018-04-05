# This python script generates the needed sub.sh, nve.sh, and job_array.sh scripts in the simulation directory i.e. /sim/
# required input: input_file
# produced output: sub.sh, nve.sh, job_array.sh

input_file='input_file'
ip = open(input_file)

molec=[]
for i, line in enumerate(ip):
    if i == 10:
        # 11th line
        num_molecs=int(line)
    if i == 12:
        molec.append(line.strip())
    if i == 13 and num_molecs > 1:
        molec.append(line.strip())
    if i == 14 and num_molecs > 2:
        molec.append(line.strip())
    if i == 15 and num_molecs > 3:
        molec.append(line.strip())
    if i == 17:
        start_config=line.strip()
    if i == 19:
        end_config=line.strip()
    if i == 21:
        sep_config=line.strip()
    if i == 23:
        timestep=line.strip()
    if i == 25: 
        num_times=line.strip()
    if i == 27: 
        nblocks=line.strip()
ip.close()

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

ja.write('SEP=%s\n' % (sep_config))
ja.write('START=%s\n' % (start_config))
ja.write('CUR=$(( MOAB_JOBARRAYINDEX*SEP - SEP + START ))\n')
ja.write('cd $PBS_O_WORKDIR\n')
ja.write('cd FILES/$CUR\n\n\n')

ja.write('module load lammps/11Aug17\n\n')

ja.write('echo Time is `date` > array_$MOAB_JOBARRAYINDEX.o\n')
ja.write('echo Directory is `pwd` >> array_$MOAB_JOBARRAYINDEX.o\n\n\n')
ja.write('mpirun lmp_mpi < in.nve -screen none\n\n\n')
for i in range(0,num_molecs):
    ja.write('python set_msd_calcs.py -inp $CUR -mol %s -ntimes %s -stp %s\n' % (molec[i], num_times, timestep))
    ja.write('./msd_rot_calc < msd_rot_calc.in\n\n')
    ja.write('python grab_press.py\n')
    ja.write('./visc_calc\n\n')

ja.write('echo Ending Time is `date` >> array_$MOAB_JOBARRAYINDEX.o\n')

# Generate NVE Input File
nve_file="nve.sh"
nve=open(nve_file, 'w')

nve.write('#MSUB -N direct_calc_nve\n')
nve.write('#MSUB -q thompson\n')
nve.write('#MSUB -j oe\n')
nve.write('#MSUB -d ./\n')
nve.write('#MSUB -l nodes=1:ppn=10:intel,mem=5gb,walltime=6:00:00\n\n\n')

nve.write('SEP=%s\n' % (sep_config))
nve.write('START=%s\n' % (start_config))
nve.write('CUR=$(( MOAB_JOBARRAYINDEX*SEP - SEP + START ))\n')
nve.write('cd $PBS_O_WORKDIR\n')
nve.write('cd FILES/$CUR\n\n\n')

nve.write('module load lammps/11Aug17\n\n')

nve.write('echo Time is `date` > array_$MOAB_JOBARRAYINDEX.o\n')
nve.write('echo Directory is `pwd` >> array_$MOAB_JOBARRAYINDEX.o\n\n\n')
nve.write('mpirun lmp_mpi < in.nve -screen none\n\n\n')
for i in range(0,num_molecs):
    nve.write('python set_msd_calcs.py -inp ${PWD##*/} -mol %s -ntimes %s -stp %s\n' % (molec[i], num_times, timestep))
    nve.write('./msd_rot_calc < msd_rot_calc.in\n\n')

nve.write('python grab_press.py\n\n')
nve.write('./visc_calc\n')

nve.write('echo Ending Time is `date` >> array_$MOAB_JOBARRAYINDEX.o\n')


nve.close()
