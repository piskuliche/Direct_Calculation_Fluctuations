#MSUB -N direct_calc_nve
#MSUB -q sixhour
#MSUB -j oe
#MSUB -d ./
#MSUB -l nodes=1:ppn=10:ib,mem=5gb,walltime=6:00:00
#MSUB -t AAA-BBB

SEP=1000
START=1001000
CUR=$(( MOAB_JOBARRAYINDEX*SEP - SEP + START ))
cd $PBS_O_WORKDIR
cd $CUR


module load legacy
module load intel_compiler_2016
module load intel_mpi_intel64/5.1.2.150
module load lammps/20151207

echo Time is `date` > array_$MOAB_JOBARRAYINDEX.o
echo Directory is `pwd` >> array_$MOAB_JOBARRAYINDEX.o

mpirun lmp_g++ < in.nve -screen none

echo Ending Time is `date` >> array_$MOAB_JOBARRAYINDEX.o

