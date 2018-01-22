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
cd FILES/$CUR


module load legacy
module load intel_compiler_2016
module load intel_mpi_intel64/5.1.2.150
module load lammps/11Aug17

echo Time is `date` > array_$MOAB_JOBARRAYINDEX.o
echo Directory is `pwd` >> array_$MOAB_JOBARRAYINDEX.o
mpirun lmp_mpi < in.nve -screen none
python set_msd_calcs.py -inp $CUR -mol acn
./msd_rot_calc < msd_rot_calc.in
python set_msd_calcs.py -inp $CUR -mol co2
./msd_rot_calc < msd_rot_calc.in


echo Ending Time is `date` >> array_$MOAB_JOBARRAYINDEX.o

